/*jslint browser:true devel:true plusplus:true indent:4 continue:true nomen:true*/
/*global GLmol, THREE */

'use strict';
/*  ProteinSurface.js by biochem_fan

Ported and modified for Javascript based on EDTSurf,
  whose license is as follows.

Permission to use, copy, modify, and distribute this program for any
purpose, with or without fee, is hereby granted, provided that this
copyright notice and the reference information appear in all copies or
substantial portions of the Software. It is provided "as is" without
express or implied warranty. 

Reference:
http://zhanglab.ccmb.med.umich.edu/EDTSurf/
D. Xu, Y. Zhang (2009) Generating Triangulated Macromolecular Surfaces
by Euclidean Distance Transform. PLoS ONE 4(12): e8140.

========

TODO: Improved performance on Firefox
      Reduce memory consumption
      Refactor!

      - think about using THREE.MarchingCubes
*/



// webworker hacks
var GLmol = GLmol || function () {};
if (!this.console) {
  var console = {
    log: function(){}
  };
}

if (typeof importScripts !== 'undefined') {
   var window = this;
   importScripts("three.min.js");
}


GLmol.prototype.generateMesh = function (group, atomlist, type, wireframe, wireframeLinewidth, async) {
    wireframe = wireframe || false;
    var atomsToShow, extent, expandedExtent, extendedAtoms, ps, mat, mesh;
    wireframeLinewidth = wireframeLinewidth || 1;

    var proteinSurfaceDone = function (surfaceGeo) {
        this.surfaceGeo = surfaceGeo;
        mesh = this.getLambertMesh(this.surfaceGeo, {
            vertexColors: THREE.VertexColors,
            wireframe: wireframe,
            wireframeLinewidth: wireframeLinewidth,
            opacity: 0.8,
            transparent: true
        })
        mesh.doubleSided = true;
        group.add(mesh);
        this.show();
    }.bind(this);

    if (!this.surfaceGeo || this.meshType !== type) {
        atomsToShow = this.removeSolvents(atomlist);
        extent = this.getExtent(atomsToShow);
        expandedExtent = [[extent[0][0] - 4, extent[0][1] - 4, extent[0][2] - 4],
                          [extent[1][0] + 4, extent[1][1] + 4, extent[1][2] + 4]];
        extendedAtoms = this.removeSolvents(this.getAtomsWithin(this.atoms, expandedExtent));
        this.meshType = type;


        if (async) {
            var worker = new SharedWorker("js/ProteinSurface4.js");
            worker.port.onmessage = function (event) {
                console.log("done");
                proteinSurfaceDone(event.data);
            }
            worker.port.start()

            worker.postMessage = worker.webkitPostMessage || worker.postMessage;

            worker.port.postMessage([expandedExtent, type, this.atoms, extendedAtoms, atomsToShow]);
        } else {
            var ps = new ProteinSurface(expandedExtent, type, this.atoms, extendedAtoms);
            proteinSurfaceDone(ps.getModel(this.atoms, atomsToShow));
        }
    }

};


this.addEventListener('connect', function(e) {
    var port = e.ports[0];  
    port.addEventListener('message', function(e) {
        var expandedExtent = e.data[0],
            type = e.data[1],
            atoms = e.data[2],
            extendedAtoms = e.data[3],
            atomsToShow = e.data[4];
        var ps = new ProteinSurface(expandedExtent, type, atoms, extendedAtoms);
        port.postMessage( ps.getModel(atoms, atomsToShow) )
    });
});


function ProteinSurface(expandedExtent, type, atoms, extendedAtoms) {
    console.profile("constructor")
    this.initparm(expandedExtent, (type === ProteinSurface.VDW) ? false : true);
    this.fillvoxels(atoms, extendedAtoms);
    this.buildboundary();
    if (type === ProteinSurface.MS || type === ProteinSurface.SES) {
        this.fastdistancemap();
    }
    if (type === ProteinSurface.SES) {
        this.boundingatom(false);
        this.fillvoxelswaals(atoms, extendedAtoms);
    }
    this.marchingcube(type);
    this.laplaciansmooth(1);
    console.profileEnd("constructor")
};

ProteinSurface.VDW = 1;
ProteinSurface.SES = 2;
ProteinSurface.SAS = 3;
ProteinSurface.MS = 4;

ProteinSurface.prototype.inarray = [];
ProteinSurface.prototype.outarray = [];
ProteinSurface.prototype.boxLength = 128;
ProteinSurface.prototype.probeRadius = 1.4;
ProteinSurface.prototype.scaleFactor = 1;
ProteinSurface.prototype.rasrad = [1.90, 1.88, 1.63, 1.48, 1.78, 1.2, 1.87, 1.96, 1.63, 0.74, 1.8, 1.48, 1.2]; //liang
//             Calpha   c    n    o    s   h   p   Cbeta  ne  fe  other ox  hx

ProteinSurface.prototype.depty = [] // new Array(13);
ProteinSurface.prototype.widxz = [] //new Int32Array(13);
ProteinSurface.prototype.fixsf = 2;
ProteinSurface.prototype.nb = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1],
           [1, 1, 0], [1, -1, 0], [-1, 1, 0], [-1, -1, 0], [1, 0, 1], [1, 0, -1],
           [-1, 0, 1], [-1, 0, -1], [0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1],
           [1, 1, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1], [1, -1, -1], [-1, -1, 1], [-1, 1, -1], [-1, -1, -1]];

ProteinSurface.prototype.getModel = function (atoms, atomlist) {
    console.profile("getModel");
    var i,
        lim,
        v = [],
        vertices = this.verts,
        atomsToShow = {},
        scaleFactor = this.scaleFactor,
        ptranx = this.ptranx,
        ptrany = this.ptrany,
        ptranz = this.ptranz,
        geo = new THREE.Geometry(),
        f,
        a,
        b,
        c;
    console.log("all vertices & faces", this.verts.length, this.faces.length);

    atomsToShow = {};

    for (i = 0, lim = atomlist.length; i < lim; i++) {
        atomsToShow[i] = true;
    }
    for (i = 0; i < this.verts.length; i++) {
        vertices[i].x = vertices[i].x / scaleFactor - ptranx;
        vertices[i].y = vertices[i].y / scaleFactor - ptrany;
        vertices[i].z = vertices[i].z / scaleFactor - ptranz;
        v.push(vertices[i]);
    }

    geo.vertices = v;
    for (i = 0; i < this.faces.length; i++) {
        f = this.faces[i];
        a = vertices[f.a].atomid;
        b = vertices[f.b].atomid;
        c = vertices[f.c].atomid;
        if (!atomsToShow[a] && !atomsToShow[b] && !atomsToShow[c]) {
            continue;
        }
        f.vertexColors = [new THREE.Color(atoms[a].color),
                          new THREE.Color(atoms[b].color),
                          new THREE.Color(atoms[c].color)];
        geo.faces.push(f);
    }
    geo.computeFaceNormals();
    geo.computeVertexNormals();
    console.profileEnd("getModel");
    return geo;
};

ProteinSurface.prototype.laplaciansmooth = function (numiter) {
    var tps = new Array(this.verts.length),
        i,
        j,
        k,
        vertdeg = new Array(20),
        faces = this.faces,
        verts = this.verts,
        flagvert,
        wt = 1.00,
        wt2 = 0.50,
        ssign,
        outwt = 0.75 / (this.scaleFactor + 3.5); //area-preserving

    for (i = 0; i < this.verts.length; i++) {
        tps[i] = {x: 0, y: 0, z: 0};
    }
    for (i = 0; i < 20; i++) {
        vertdeg[i] = new Array(verts.length);
    }

    for (i = 0; i < this.verts.length; i++) {
        vertdeg[0][i] = 0;
    }
    for (i = 0; i < this.faces.length; i++) {
        flagvert = true;
        for (j = 0; j < vertdeg[0][faces[i].a]; j++) {
            if (faces[i].b === vertdeg[j + 1][faces[i].a]) {
                flagvert = false;
                break;
            }
        }
        if (flagvert) {
            vertdeg[0][faces[i].a]++;
            vertdeg[vertdeg[0][faces[i].a]][faces[i].a] = faces[i].b;
        }
        flagvert = true;
        for (j = 0; j < vertdeg[0][faces[i].a]; j++) {
            if (faces[i].c === vertdeg[j + 1][faces[i].a]) {
                flagvert = false;
                break;
            }
        }
        if (flagvert) {
            vertdeg[0][faces[i].a]++;
            vertdeg[vertdeg[0][faces[i].a]][faces[i].a] = faces[i].c;
        }
        //b
        flagvert = true;
        for (j = 0; j < vertdeg[0][faces[i].b]; j++) {
            if (faces[i].a === vertdeg[j + 1][faces[i].b]) {
                flagvert = false;
                break;
            }
        }
        if (flagvert) {
            vertdeg[0][faces[i].b]++;
            vertdeg[vertdeg[0][faces[i].b]][faces[i].b] = faces[i].a;
        }

        flagvert = true;

        for (j = 0; j < vertdeg[0][faces[i].b]; j++) {
            if (faces[i].c === vertdeg[j + 1][faces[i].b]) {
                flagvert = false;
                break;
            }
        }

        if (flagvert) {
            vertdeg[0][faces[i].b]++;
            vertdeg[vertdeg[0][faces[i].b]][faces[i].b] = faces[i].c;
        }
         //c
        flagvert = true;
        for (j = 0; j < vertdeg[0][faces[i].c]; j++) {
            if (faces[i].a === vertdeg[j + 1][faces[i].c]) {
                flagvert = false;
                break;
            }
        }
        if (flagvert) {
            vertdeg[0][faces[i].c]++;
            vertdeg[vertdeg[0][faces[i].c]][faces[i].c] = faces[i].a;
        }
        flagvert = true;
        for (j = 0; j < vertdeg[0][faces[i].c]; j++) {
            if (faces[i].b === vertdeg[j + 1][faces[i].c]) {
                flagvert = false;
                break;
            }
        }
        if (flagvert) {
            vertdeg[0][faces[i].c]++;
            vertdeg[vertdeg[0][faces[i].c]][faces[i].c] = faces[i].b;
        }
    }

    for (k = 0; k < numiter; k++) {
        for (i = 0; i < verts.length; i++) {
            if (vertdeg[0][i] < 3) {
                tps[i].x = verts[i].x;
                tps[i].y = verts[i].y;
                tps[i].z = verts[i].z;
            } else if (vertdeg[0][i] === 3 || vertdeg[0][i] === 4) {
                tps[i].x = 0;
                tps[i].y = 0;
                tps[i].z = 0;
                for (j = 0; j < vertdeg[0][i]; j++) {
                    tps[i].x += verts[vertdeg[j + 1][i]].x;
                    tps[i].y += verts[vertdeg[j + 1][i]].y;
                    tps[i].z += verts[vertdeg[j + 1][i]].z;
                }
                tps[i].x += wt2 * verts[i].x;
                tps[i].y += wt2 * verts[i].y;
                tps[i].z += wt2 * verts[i].z;
                tps[i].x /= wt2 + vertdeg[0][i];
                tps[i].y /= wt2 + vertdeg[0][i];
                tps[i].z /= wt2 + vertdeg[0][i];
            } else {
                tps[i].x = 0;
                tps[i].y = 0;
                tps[i].z = 0;
                for (j = 0; j < vertdeg[0][i]; j++) {
                    tps[i].x += verts[vertdeg[j + 1][i]].x;
                    tps[i].y += verts[vertdeg[j + 1][i]].y;
                    tps[i].z += verts[vertdeg[j + 1][i]].z;
                }
                tps[i].x += wt * verts[i].x;
                tps[i].y += wt * verts[i].y;
                tps[i].z += wt * verts[i].z;
                tps[i].x /= wt + vertdeg[0][i];
                tps[i].y /= wt + vertdeg[0][i];
                tps[i].z /= wt + vertdeg[0][i];
            }
        }
        for (i = 0; i < verts.length; i++) {
            verts[i].x = tps[i].x;
            verts[i].y = tps[i].y;
            verts[i].z = tps[i].z;
        }
       /*	computenorm();
        for (var i = 0; i < verts.length; i++) {
            if (verts[i].inout) ssign = 1;
            else ssign = -1;
            verts[i].x += ssign * outwt * verts[i].pn.x;
            verts[i].y += ssign * outwt * verts[i].pn.y;
            verts[i].z += ssign * outwt * verts[i].pn.z;
        }*/
    }
};

ProteinSurface.prototype.initparm = function (extent, btype) {
    var margin = 2.5,
        ptranx,
        ptrany,
        ptranz,
        threshbox,
        sfthresh,
        cutRadius;

    this.pminx = extent[0][0];
    this.pmaxx = extent[1][0];
    this.pminy = extent[0][1];
    this.pmaxy = extent[1][1];
    this.pminz = extent[0][2];
    this.pmaxz = extent[1][2];
    if (btype) {
        this.pminx -= margin;
        this.pminy -= margin;
        this.pminz -= margin;

        this.pmaxx += margin;
        this.pmaxy += margin;
        this.pmaxz += margin;

    } else {
        this.pminx -= this.probeRadius + margin;
        this.pminy -= this.probeRadius + margin;
        this.pminz -= this.probeRadius + margin;
        this.pmaxx += this.probeRadius + margin;
        this.pmaxy += this.probeRadius + margin;
        this.pmaxz += this.probeRadius + margin;
    }

    this.ptranx = -this.pminx;
    this.ptrany = -this.pminy;
    this.ptranz = -this.pminz;
    this.scaleFactor = this.pmaxx - this.pminx;
    if ((this.pmaxy - this.pminy) > this.scaleFactor) {
        this.scaleFactor = this.pmaxy - this.pminy;
    }
    if ((this.pmaxz - this.pminz) > this.scaleFactor) {
        this.scaleFactor = this.pmaxz - this.pminz;
    }
    this.scaleFactor = (this.boxLength - 1.0) / this.scaleFactor;

    this.boxLength = Math.floor(this.boxLength * this.fixsf / this.scaleFactor);
    this.scaleFactor =  this.fixsf;
    threshbox = 180; // maximum possible boxsize // FIXME magic number
    if (this.boxLength > threshbox) {
        sfthresh = threshbox / this.boxLength;
        this.boxLength = Math.floor(threshbox);
        this.scaleFactor = this.scaleFactor * sfthresh;
    }

    this.pLength = Math.ceil(this.scaleFactor * (this.pmaxx - this.pminx)) + 1;
    this.pWidth = Math.ceil(this.scaleFactor * (this.pmaxy - this.pminy)) + 1;
    this.pHeight = Math.ceil(this.scaleFactor * (this.pmaxz - this.pminz)) + 1;
    if (this.pLength > this.boxLength) {
        this.pLength = this.boxLength;
    }
    if (this.pWidth > this.boxLength) {
        this.pWidth = this.boxLength;
    }
    if (this.pHeight > this.boxLength) {
        this.pHeight = this.boxLength;
    }
    this.boundingatom(btype);
    this.cutRadius = this.probeRadius * this.scaleFactor;

    this.vp = new Array(this.pLength * this.pWidth * this.pHeight);
    console.log("Box size: ", this.pLength, this.pWidth, this.pHeight, this.vp.length);
};

ProteinSurface.prototype.boundingatom = function (btype) {
    var tradius = new Array(13),
        txz,
        tdept,
        sradius,
        idx,
        flagradius = btype,
        i,
        j,
        k,
        scaleFactor = this.scaleFactor,
        probeRadius = this.probeRadius,
        rasrad = this.rasrad,
        indx;

    for (i = 0; i < 13; i++) {
        if (!btype) {
            tradius[i] = rasrad[i] * scaleFactor + 0.5;
        } else {
            tradius[i] = (rasrad[i] + probeRadius) * scaleFactor + 0.5;
        }

        sradius = tradius[i] * tradius[i];
        this.widxz[i] = Math.floor(tradius[i]) + 1;
        this.depty[i] = new Array(this.widxz[i] * this.widxz[i]);
        indx = 0;
        for (j = 0; j < this.widxz[i]; j++) {
            for (k = 0; k < this.widxz[i]; k++) {
                txz = j * j + k * k;
                if (txz > sradius) {
                    this.depty[i][indx] = -1; // outside
                } else {
                    tdept = Math.sqrt(sradius - txz);
                    this.depty[i][indx] = Math.floor(tdept);
                }
                indx++;
            }
        }
    }
};

ProteinSurface.prototype.fillvoxels = function (atoms, atomlist) { //(int seqinit,int seqterm,bool atomtype,atom* proseq,bool bcolor)
    var i,
        vp = this.vp,
        lim,
        atom;
    for (i = 0, lim = vp.length; i < lim; i++) {
        vp[i] = {inout: false, isdone: false, isbound: false, distance: -1, atomid: -1};
    }

    for (i in atomlist) {
        if (atomlist.hasOwnProperty(i)) {
            atom = atomlist[i];
            if (!atom || atom.hetflag) { continue; }
            this.fillAtom(atom, atoms);
        }
    }

    for (i = 0, lim = vp.length; i < lim; i++) {
        if (vp[i].inout) {
            vp[i].isdone = true;
        }
    }

    this.vp = vp;
    for (i = 0, lim = vp.length; i < lim; i++) {
        if (vp[i].inout) {
            vp[i].isdone = true;
        }
    }
};

ProteinSurface.prototype.getAtomType = function (atom) {
  //atom = {
    //'CA': 0,
    //'C': 1,
    //'0', 3,
    //'N': 2,
    //'FE': 9,
    //'H': 5
  //}
  //elem = {
    //'C': 1,
    //'O': 11,
    //'N': 8,
    //'S': 4,
    //'P': 6,
    //'H': 1,
  //}
    var at = 10;
    if (atom.atom === 'CA') {
        at = 0;
    } else if (atom.atom === 'C') {
        at = 1;
    } else if (atom.elem === 'C') {
        at = 7;
    } else if (atom.atom === '0') {
        at = 3;
    } else if (atom.elem === 'O') {
        at = 11;
    } else if (atom.atom === 'N') {
        at = 2;
    } else if (atom.elem === 'N') {
        at = 8;
    } else if (atom.elem === 'S') {
        at = 4;
    } else if (atom.elem === 'P') {
        at = 6;
    } else if (atom.atom === 'FE') {
        at = 9;
    } else if (atom.atom === 'H') {
        at = 5;
    } else if (atom.elem === 'H') {
        at = 1;
    }
    return at;
};

ProteinSurface.prototype.fillAtom = function (atom, atoms) {
    var cx,
        cy,
        cz,
        ox,
        oy,
        oz,
        scaleFactor = this.scaleFactor,
        ptranx = this.ptranx,
        ptrany = this.ptranz,
        ptranz = this.ptranz,
        at = this.getAtomType(atom),
        widxz = this.widxz,
        depty = this.depty,
        pHeight = this.pHeight,
        pWidth = this.pWidth,
        pLength = this.pLength,
        nind = 0,
        cnt = 0,
        i,
        j,
        ii,
        jj,
        kk,
        mi,
        mk,
        k,
        mj,
        si,
        sj,
        sk,
        vpSISJSK,
        atom2;


    cx = Math.floor(0.5 + scaleFactor * (atom.x + ptranx));
    cy = Math.floor(0.5 + scaleFactor * (atom.y + ptrany));
    cz = Math.floor(0.5 + scaleFactor * (atom.z + ptranz));


    for (i = 0; i < widxz[at]; i++) {
        for (j = 0; j < widxz[at]; j++) {
            if (depty[at][nind] !== -1) {
                for (ii = -1; ii < 2; ii++) {
                    for (jj = -1; jj < 2; jj++) {
                        for (kk = -1; kk < 2; kk++) {
                            if (ii !== 0 && jj !== 0 && kk !== 0) {
                                mi = ii * i;
                                mk = kk * j;
                                for (k = 0; k <= depty[at][nind]; k++) {
                                    mj = k * jj;
                                    si = cx + mi;
                                    sj = cy + mj;
                                    sk = cz + mk;
                                    if (si < 0 || sj < 0 || sk < 0 || si >= pLength || sj >= pWidth || sk >= pHeight) {
                                        continue;
                                    }
                                    vpSISJSK = this.vp[si * pWidth * pHeight + sj * pHeight + sk];
                                    if (false) { // !bcolor
                                        vpSISJSK.inout = true;
                                    } else { // color
                                        if (!vpSISJSK.inout) {
                                            vpSISJSK.inout = true;
                                            vpSISJSK.atomid = atom.serial;
                                        } else if (vpSISJSK.inout) {
                                            atom2 = atoms[vpSISJSK.atomid];
                                            ox = Math.floor(0.5 + scaleFactor * (atom2.x + ptranx));
                                            oy = Math.floor(0.5 + scaleFactor * (atom2.y + ptrany));
                                            oz = Math.floor(0.5 + scaleFactor * (atom2.z + ptranz));
                                            if (mi * mi + mj * mj + mk * mk < ox * ox + oy * oy + oz * oz) {
                                                vpSISJSK.atomid = atom.serial;
                                            }
                                        }
                                    }
                                }//k
                            }//if
                        }//kk	
                    }//jj
                }//ii
            }//if
            nind++;
        }//j
    }//i
};

ProteinSurface.prototype.fillvoxelswaals = function (atoms, atomlist) {
    var i,
        vp = this.vp,
        lim = vp.length,
        atom;

    for (i = 0, lim; i < lim; i++) {
        vp[i].isdone = false;
    }

    for (i in atomlist) {
        if (atomlist.hasOwnProperty(i)) {
            atom = atomlist[i];
            if (!atom || atom.hetflag) {
                continue;
            }

            this.fillAtomWaals(atom, atoms);
        }
    }
};

ProteinSurface.prototype.fillAtomWaals = function (atom, atoms) {
    var cx,
        cy,
        cz,
        ox,
        oy,
        oz,
        nind = 0,
        scaleFactor = this.scaleFactor,
        ptranx = this.ptranx,
        ptrany = this.ptranz,
        ptranz = this.ptranz,
        at = this.getAtomType(atom),
        widxz = this.widxz,
        depty = this.depty,
        pHeight = this.pHeight,
        pWidth = this.pWidth,
        pLength = this.pLength,
        vp = this.vp,
        i,
        j,
        k,
        ii,
        jj,
        kk,
        mi,
        mk,
        mj,
        si,
        sj,
        sk,
        vpSISJSK,
        atom2;

    cx = Math.floor(0.5 + scaleFactor * (atom.x + ptranx));
    cy = Math.floor(0.5 + scaleFactor * (atom.y + ptrany));
    cz = Math.floor(0.5 + scaleFactor * (atom.z + ptranz));


    for (i = 0; i < widxz[at]; i++) {
        for (j = 0; j < widxz[at]; j++) {
            if (depty[at][nind] !== -1) {
                for (ii = -1; ii < 2; ii++) {
                    for (jj = -1; jj < 2; jj++) {
                        for (kk = -1; kk < 2; kk++) {
                            if (ii !== 0 && jj !== 0 && kk !== 0) {
                                mi = ii * i;
                                mk = kk * j;
                                for (k = 0; k <= depty[at][nind]; k++) {
                                    mj = k * jj;
                                    si = cx + mi;
                                    sj = cy + mj;
                                    sk = cz + mk;
                                    if (si < 0 || sj < 0 || sk < 0) { continue; }
                                    vpSISJSK = vp[si * pWidth * pHeight + sj * pHeight + sk];
                                    if (false) { //(!bcolor) FIXME
                                        vpSISJSK.isdone = true;
                                        continue;
                                    } else {
                                        if (vpSISJSK.isdone === false) {
                                            vpSISJSK.isdone = true;
                                            vpSISJSK.atomid = atom.serial;
                                        } else if (vpSISJSK.isdone) {
                                            atom2 = atoms[vpSISJSK.atomid];
                                            ox = Math.floor(0.5 + scaleFactor * (atom2.x + ptranx));
                                            oy = Math.floor(0.5 + scaleFactor * (atom2.y + ptrany));
                                            oz = Math.floor(0.5 + scaleFactor * (atom2.z + ptranz));
                                            if (mi * mi + mj * mj + mk * mk < ox * ox + oy * oy + oz * oz) {
                                                vpSISJSK.atomid = atom.serial;
                                            }
                                        }
                                    }//else
                                 //k
                                }//if
                            }//kk	
                        }//jj
                    }//ii
                }//if
                nind++;
            }//j
        }//i
    }
};


ProteinSurface.prototype.buildboundary = function () {
    var vp = this.vp,
        pHeight = this.pHeight,
        pWidth = this.pWidth,
        pLength = this.pLength,
        i,
        j,
        k,
        ii,
        vpIJK,
        flagbound,
        ti,
        tj,
        tk;

    for (i = 0; i < pLength; i++) {
        for (j = 0; j < pHeight; j++) {
            for (k = 0; k < pWidth; k++) {
                vpIJK = vp[i * pWidth * pHeight + k * pHeight + j];
                if (vpIJK.inout) {
                    flagbound = false;
                    ii = 0;
                    while (!flagbound && ii < 26) {
                        ti = i + this.nb[ii][0];
                        tj = j + this.nb[ii][2];
                        tk = k + this.nb[ii][1];
                        if (ti > -1 && ti < pLength
                                && tk > -1 && tk < pWidth
                                && tj > -1 && tj < pHeight
                                && !vp[ti * pWidth * pHeight + tk * pHeight + tj].inout) {
                            vpIJK.isbound = true;
                            flagbound = true;
                        } else {
                            ii++;
                        }
                    }
                }
            }
        }
    }
};

ProteinSurface.prototype.fastdistancemap = function () {
    var positin,
        vp = this.vp,
        pHeight = this.pHeight,
        pWidth = this.pWidth,
        pLength = this.pLength,
        positout,
        eliminate,
        certificate,
        totalsurfacevox = 0,
        totalinnervox = 0,
        boundPoint = new Array(pLength),
        i,
        a,
        j,
        b,
        k,
        vpIJK,
        vptmp,
        cutRadis = this.cutRadis,
        scaleFactor = this.scaleFactor,
        cutsf = scaleFactor - 0.5;

    for (i = 0; i < pLength; i++) {
        a = new Array(pWidth);
        for (j = 0; j < pWidth; j++) {
            b = new Array(pHeight);
            for (k = 0; k < pHeight; k++) {
                b[k] = {ix: 0, iy: 0, iz: 0};
            }
            a[j] = b;
        }
        boundPoint[i] = a;
    }

    for (i = 0; i < pLength; i++) {
        for (j = 0; j < pWidth; j++) {
            for (k = 0; k < pHeight; k++) {
                vpIJK = vp[i * pWidth * pHeight + j * pHeight + k];
                vpIJK.isdone = false;
                if (vpIJK.inout) {
                    if (vpIJK.isbound) {
                        totalsurfacevox++;
                        boundPoint[i][j][k].ix = i;
                        boundPoint[i][j][k].iy = j;
                        boundPoint[i][j][k].iz = k;
                        vpIJK.distance = 0;
                        vpIJK.isdone = true;
                    } else {
                        totalinnervox++;
                    }
                }
            }
        }
    }

    positin = 0;
    positout = 0;

    for (i = 0; i < pLength; i++) {
        for (j = 0; j < pWidth; j++) {
            for (k = 0; k < pHeight; k++) {
                vpIJK = vp[i * pWidth * pHeight + j * pHeight + k];
                if (vpIJK.isbound) {
                    this.inarray.push({ix: i, iy: j, iz: k});
                    positin++;
                    vpIJK.isbound = false;
                }
            }
        }
    }

    do {
        positout = this.fastoneshell(positin, boundPoint);
        positin = 0;
        this.inarray = [];
        for (i = 0; i < positout; i++) {
            vptmp = vp[pWidth * pHeight * this.outarray[i].ix + pHeight * this.outarray[i].iy + this.outarray[i].iz];
            vptmp.isbound = false;
            if (vptmp.distance <= 1.02 * cutRadis) {
                this.inarray.push({ix: outarray[i].ix, iy: outarray[i].iy, iz: outarray[i].iz});
                //            inarray[positin].ix=outarray[i].ix;
                //            inarray[positin].iy=outarray[i].iy;
                //            inarray[positin].iz=outarray[i].iz;
                positin++;
            }
        }
    } while (positin !== 0);

    if (cutsf < 0) {
        cutsf = 0;
    }
    for (i = 0; i < pLength; i++) {
        for (j = 0; j < pWidth; j++) {
            for (k = 0; k < pHeight; k++) {
                vpIJK = vp[i * pWidth * pHeight + j * pHeight + k];
                vpIJK.isbound = false;
                //ses solid
                if (vpIJK.inout) {
                    if (!vpIJK.isdone || (vpIJK.isdone && vpIJK.distance >= cutRadis - 0.50 / (0.1 + cutsf))) {
                        vpIJK.isbound = true;
                        //new add
                        //                  if (vpIJK.isdone)
                        //                  VPIJK.atomid=vp[boundPoint[i][j][k].ix][boundPoint[i][j][k].iy][boundPoint[i][j][k].iz].atomid;
                    }
                }
            }
        }
    }
    this.inarray = [];
    this.outarray = [];
};

ProteinSurface.prototype.fastoneshell = function (number, boundPoint) { //(int* innum,int *allocout,voxel2 ***boundPoint, int* outnum, int *elimi)
    var positout = 0,
        tx,
        ty,
        tz,
        dx,
        dy,
        dz,
        square,
        tnv = {ix: -1, iy: -1, iz: -1},
        i,
        j,
        nb = this.nb,
        vpTNV,
        pHeight = this.pHeight,
        pWidth = this.pWidth,
        pLength = this.pLength,
        vp = this.vp;

    if (number === 0) {
        return 0;
    }

    for (i = 0; i < number; i++) {
        tx = this.inarray[i].ix;
        ty = this.inarray[i].iy;
        tz = this.inarray[i].iz;

        for (j = 0; j < 6; j++) {
            tnv.ix = tx + nb[j][0];
            tnv.iy = ty + nb[j][1];
            tnv.iz = tz + nb[j][2];
            vpTNV = vp[tnv.ix * pWidth * pHeight + pHeight * tnv.iy + tnv.iz];
            if (tnv.ix < pLength && tnv.ix > -1 &&
                    tnv.iy < pWidth && tnv.iy > -1 &&
                    tnv.iz < pHeight && tnv.iz > -1 &&
                    vpTNV.inout && !vpTNV.isdone) {
                boundPoint[tnv.ix][tnv.iy][tz + nb[j][2]].ix = boundPoint[tx][ty][tz].ix;
                boundPoint[tnv.ix][tnv.iy][tz + nb[j][2]].iy = boundPoint[tx][ty][tz].iy;
                boundPoint[tnv.ix][tnv.iy][tz + nb[j][2]].iz = boundPoint[tx][ty][tz].iz;
                dx = tnv.ix - boundPoint[tx][ty][tz].ix;
                dy = tnv.iy - boundPoint[tx][ty][tz].iy;
                dz = tnv.iz - boundPoint[tx][ty][tz].iz;
                square = dx * dx + dy * dy + dz * dz;
                vpTNV.distance = Math.sqrt(square);
                vpTNV.isdone = true;
                vpTNV.isbound = true;
                this.outarray.push({ix: tnv.ix, iy: tnv.iy, iz: tnv.iz});
                positout++;
            } else if (tnv.ix < pLength && tnv.ix > -1 &&
                    tnv.iy < pWidth && tnv.iy > -1 &&
                    tnv.iz < pHeight && tnv.iz > -1 &&
                       vpTNV.inout && vpTNV.isdone) {
                dx = tnv.ix - boundPoint[tx][ty][tz].ix;
                dy = tnv.iy - boundPoint[tx][ty][tz].iy;
                dz = tnv.iz - boundPoint[tx][ty][tz].iz;
                square = dx * dx + dy * dy + dz * dz;
                square = Math.sqrt(square);
                if (square < vpTNV.distance) {
                    boundPoint[tnv.ix][tnv.iy][tnv.iz].ix = boundPoint[tx][ty][tz].ix;
                    boundPoint[tnv.ix][tnv.iy][tnv.iz].iy = boundPoint[tx][ty][tz].iy;
                    boundPoint[tnv.ix][tnv.iy][tnv.iz].iz = boundPoint[tx][ty][tz].iz;
                    vpTNV.distance = square;
                    if (!vpTNV.isbound) {
                        vpTNV.isbound = true;
                        this.outarray.push({ix: tnv.ix, iy: tnv.iy, iz: tnv.iz});
                        positout++;
                    }
                }
            }
        }
    }

    console.log("part1", positout);

    for (i = 0; i < number; i++) {
        tx = this.inarray[i].ix;
        ty = this.inarray[i].iy;
        tz = this.inarray[i].iz;
        for (j = 6; j < 18; j++) {
            tnv.ix = tx + nb[j][0];
            tnv.iy = ty + nb[j][1];
            tnv.iz = tz + nb[j][2];
            vpTNV = vp[tnv.ix * pWidth * pHeight + pHeight * tnv.iy + tnv.iz];

            if (tnv.x < pLength && tnv.ix > -1 &&
                    tnv.iy < pWidth && tnv.iy > -1 &&
                    tnv.iz < pHeight && tnv.iz  > -1 &&
                    vpTNV.inout && !vpTNV.isdone) {
                boundPoint[tnv.ix][tnv.iy][tz + nb[j][2]].ix = boundPoint[tx][ty][tz].ix;
                boundPoint[tnv.ix][tnv.iy][tz + nb[j][2]].iy = boundPoint[tx][ty][tz].iy;
                boundPoint[tnv.ix][tnv.iy][tz + nb[j][2]].iz = boundPoint[tx][ty][tz].iz;
                dx = tnv.ix - boundPoint[tx][ty][tz].ix;
                dy = tnv.iy - boundPoint[tx][ty][tz].iy;
                dz = tnv.iz - boundPoint[tx][ty][tz].iz;
                square = dx * dx + dy * dy + dz * dz;
                vpTNV.distance = Math.sqrt(square);
                vpTNV.isdone = true;
                vpTNV.isbound = true;
                this.outarray.push({ix: tnv.ix, iy: tnv.iy, iz: tnv.iz});
                positout++;
            } else if (tnv.ix < pLength && tnv.ix > -1 &&
                        tnv.iy < pWidth && tnv.iy > -1 &&
                        tnv.iz < pHeight && tnv.iz > -1 &&
                        vpTNV.inout && vpTNV.isdone) {
                dx = tnv.ix - boundPoint[tx][ty][tz].ix;
                dy = tnv.iy - boundPoint[tx][ty][tz].iy;
                dz = tnv.iz - boundPoint[tx][ty][tz].iz;
                square = Math.sqrt(dx * dx + dy * dy + dz * dz);
                if (square < vpTNV.distance) {
                    boundPoint[tnv.ix][tnv.iy][tnv.iz].ix = boundPoint[tx][ty][tz].ix;
                    boundPoint[tnv.ix][tnv.iy][tnv.iz].iy = boundPoint[tx][ty][tz].iy;
                    boundPoint[tnv.ix][tnv.iy][tnv.iz].iz = boundPoint[tx][ty][tz].iz;
                    vpTNV.distance = square;
                    if (!vpTNV.isbound) {
                        vpTNV.isbound = true;
                        this.outarray.push({ix: tnv.ix, iy: tnv.iy, iz: tnv.iz});
                        positout++;
                    }
                }
            }
        }
    }

    console.log("part2", positout);

    for (i = 0; i < number; i++) {
        tx = this.inarray[i].ix;
        ty = this.inarray[i].iy;
        tz = this.inarray[i].iz;
        for (j = 18; j < 26; j++) {
            tnv.ix = tx + nb[j][0];
            tnv.iy = ty + nb[j][1];
            tnv.iz = tz + nb[j][2];
            vpTNV = vp[tnv.ix * pWidth * pHeight + pHeight * tnv.iy + tnv.iz];

            if (tnv.ix < pLength && tnv.ix > -1 &&
                    tnv.iy < pWidth && tnv.iy > -1 &&
                    tnv.iz < pHeight && tnv.iz > -1 &&
                    vpTNV.inout && !vpTNV.isdone) {
                boundPoint[tnv.ix][tnv.iy][tz + nb[j][2]].ix = boundPoint[tx][ty][tz].ix;
                boundPoint[tnv.ix][tnv.iy][tz + nb[j][2]].iy = boundPoint[tx][ty][tz].iy;
                boundPoint[tnv.ix][tnv.iy][tz + nb[j][2]].iz = boundPoint[tx][ty][tz].iz;
                dx = tnv.ix - boundPoint[tx][ty][tz].ix;
                dy = tnv.iy - boundPoint[tx][ty][tz].iy;
                dz = tnv.iz - boundPoint[tx][ty][tz].iz;
                square = dx * dx + dy * dy + dz * dz;
                vpTNV.distance = Math.sqrt(square);
                vpTNV.isdone = true;
                vpTNV.isbound = true;
                this.outarray.push({ix: tnv.ix, iy: tnv.iy, iz: tnv.iz});
                positout++;
            } else if (tnv.ix < pLength && tnv.ix > -1 &&
                       tnv.iy < pWidth && tnv.iy > -1 &&
                       tnv.iz < pHeight && tnv.iz > -1 &&
                       vpTNV.inout && vpTNV.isdone) {
                dx = tnv.ix - boundPoint[tx][ty][tz].ix;
                dy = tnv.iy - boundPoint[tx][ty][tz].iy;
                dz = tnv.iz - boundPoint[tx][ty][tz].iz;
                square = Math.sqrt(dx * dx + dy * dy + dz * dz);
                if (square < vpTNV.distance) {
                    boundPoint[tnv.ix][tnv.iy][tnv.iz].ix = boundPoint[tx][ty][tz].ix;
                    boundPoint[tnv.ix][tnv.iy][tnv.iz].iy = boundPoint[tx][ty][tz].iy;
                    boundPoint[tnv.ix][tnv.iy][tnv.iz].iz = boundPoint[tx][ty][tz].iz;
                    vpTNV.distance = square;
                    if (!vpTNV.isbound) {
                        vpTNV.isbound = true;
                        this.outarray.push({ix: tnv.ix, iy: tnv.iy, iz: tnv.iz});
                        positout++;
                    }
                }
            }
        }
    }

    console.log("part3", positout);
    return positout;
};

ProteinSurface.prototype.marchingcubeinit = function (stype) {
    var i,
        vp = this.vp,
        lim;

    for (i = 0, lim = vp.length; i < lim; i++) {
        if (stype === 3) { // vdw
            vp[i].isbound = false;
        } else if (stype === 4) { // ses
            vp[i].isdone = false;
            if (vp[i].isbound) {
                vp[i].isdone = true;
            }
            vp[i].isbound = false;
        } else if (stype === 2) { // after vdw
            if (vp[i].isbound && vp[i].isdone) {
                vp[i].isbound = false;
            } else if (vp[i].isbound && !vp[i].isdone) {
                vp[i].isdone = true;
            }
        } else if (stype === 3) { //sas
            vp[i].isbound = false;
        }
    }
};

ProteinSurface.prototype.marchingcube = function (stype) {
    var pHeight = this.pHeight,
        pWidth = this.pWidth,
        pLength = this.pLength,
        vertseq = new Array(pLength),
        i,
        a = new Array(pWidth),
        j,
        b,
        k,
        verts = [],
        faces = [],
        sumtype,
        ii,
        jj,
        kk,
        tp = new Array(6),
        vp = this.vp,
        vp000,
        vp001,
        vp010,
        vp011,
        vp100,
        vp101,
        vp110,
        vp111;

    this.marchingcubeinit(stype);

    for (i = 0; i < pLength; i++) {
        a = new Array(pWidth);
        for (j = 0; j < pWidth; j++) {
            b = new Array(pHeight);
            for (k = 0; k < pHeight; k++) {
                b[k] = -1;
            }
            a[j] = b;
        }
        vertseq[i] = a;
    }
//i verts
    //(4 * (pHeight * pLength + pWidth * pLength + pHeight * pWidth)); // CHECK: Is this enough?
// //  for (var i = 0, lim = verts.length; i < lim; i++) verts[i] = new THREE.Vector3(0, 0, 0);
//   faces = [] //12 * (pHeight * pLength + pWidth * pLength + pHeight * pWidth)); // CHECK! 4
// for (var i = 0, lim = faces.length; i < lim; i++) faces[i] = new THREE.Face3(0, 0, 0);	

    for (i = 0; i < 6; i++) {
        tp[i] = new Array(3);
    }

    //face1
    for (i = 0; i < 1; i++) {
        for (j = 0; j < pWidth - 1; j++) {
            for (k = 0; k < pHeight - 1; k++) {
                vp000 = vp[pWidth * pHeight * i + pHeight * j + k].isdone;
                vp001 = vp[pWidth * pHeight * i + pHeight * j + k + 1].isdone;
                vp010 = vp[pWidth * pHeight * i + pHeight * (j + 1) + k].isdone;
                vp011 = vp[pWidth * pHeight * i + pHeight * (j + 1) + k + 1].isdone;
                vp100 = vp[pWidth * pHeight * (i + 1) + pHeight * j + k].isdone;
                vp101 = vp[pWidth * pHeight * (i + 1) + pHeight * j + k + 1].isdone;
                vp110 = vp[pWidth * pHeight * (i + 1) + pHeight * (j + 1) + k].isdone;
                vp111 = vp[pWidth * pHeight * (i + 1) + pHeight * (j + 1) + k + 1].isdone;

                if (vp000 && vp010 && vp011 && vp001) {
                    tp[0][0] = i;
                    tp[0][1] = j;
                    tp[0][2] = k;
                    tp[1][0] = i;
                    tp[1][1] = j + 1;
                    tp[1][2] = k;
                    tp[2][0] = i;
                    tp[2][1] = j + 1;
                    tp[2][2] = k + 1;
                    tp[3][0] = i;
                    tp[3][1] = j;
                    tp[3][2] = k + 1;
                    for (ii = 0; ii < 4; ii++) {
                        if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                            vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                            verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                            
                        }
                    }
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]]));
                    
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[3][0]][tp[3][1]][tp[3][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                    
                } else if ((vp000 && vp010 && vp011) || (vp010 && vp011 && vp001) || (vp011 && vp001 && vp000) || (vp001 && vp000 && vp010)) {
                    if (vp000 && vp010 && vp011) {
                        tp[0][0] = i;
                        tp[0][1] = j;
                        tp[0][2] = k;
                        tp[1][0] = i;
                        tp[1][1] = j + 1;
                        tp[1][2] = k;
                        tp[2][0] = i;
                        tp[2][1] = j + 1;
                        tp[2][2] = k + 1;
                    } else if (vp010 && vp011 && vp001) {
                        tp[0][0] = i;
                        tp[0][1] = j + 1;
                        tp[0][2] = k;
                        tp[1][0] = i;
                        tp[1][1] = j + 1;
                        tp[1][2] = k + 1;
                        tp[2][0] = i;
                        tp[2][1] = j;
                        tp[2][2] = k + 1;

                    } else if (vp011 && vp001 && vp000) {
                        tp[0][0] = i;
                        tp[0][1] = j + 1;
                        tp[0][2] = k + 1;
                        tp[1][0] = i;
                        tp[1][1] = j;
                        tp[1][2] = k + 1;
                        tp[2][0] = i;
                        tp[2][1] = j;
                        tp[2][2] = k;
                    } else if (vp001 && vp000 && vp010) {
                        tp[0][0] = i;
                        tp[0][1] = j;
                        tp[0][2] = k + 1;
                        tp[1][0] = i;
                        tp[1][1] = j;
                        tp[1][2] = k;
                        tp[2][0] = i;
                        tp[2][1] = j + 1;
                        tp[2][2] = k;

                    }
                    for (ii = 0; ii < 3; ii++) {
                        if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                            vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                            verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                            
                        }
                    }
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]]));
                    
                }
            }
        }
    }
    console.log(1);
    //face3
    for (i = 0; i < pLength - 1; i++) {
        for (j = 0; j < 1; j++) {
            for (k = 0; k < pHeight - 1; k++) {
                vp000 = vp[pWidth * pHeight * i + pHeight * j + k].isdone;
                vp001 = vp[pWidth * pHeight * i + pHeight * j + k + 1].isdone;
                vp100 = vp[pWidth * pHeight * (i + 1) + pHeight * j + k].isdone;
                vp101 = vp[pWidth * pHeight * (i + 1) + pHeight * j + k + 1].isdone;

                if (vp000 && vp100 && vp101 && vp001) {
                    tp[0][0] = i;
                    tp[0][1] = j;
                    tp[0][2] = k;
                    tp[1][0] = i + 1;
                    tp[1][1] = j;
                    tp[1][2] = k;
                    tp[2][0] = i + 1;
                    tp[2][1] = j;
                    tp[2][2] = k + 1;
                    tp[3][0] = i;
                    tp[3][1] = j;
                    tp[3][2] = k + 1;
                    for (ii = 0; ii < 4; ii++) {
                        if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                            vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                            verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                            
                        }
                    }
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]],  vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                    
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]],  vertseq[tp[3][0]][tp[3][1]][tp[3][2]]));
                    
                } else if ((vp000 && vp100 && vp101) || (vp100 && vp101 && vp001) || (vp101 && vp001 && vp000) || (vp001 && vp000 && vp100)) {
                    if (vp000 && vp100 && vp101) {
                        tp[0][0] = i;
                        tp[0][1] = j;
                        tp[0][2] = k;
                        tp[1][0] = i + 1;
                        tp[1][1] = j;
                        tp[1][2] = k;
                        tp[2][0] = i + 1;
                        tp[2][1] = j;
                        tp[2][2] = k + 1;

                    } else if (vp100 && vp101 && vp001) {
                        tp[0][0] = i + 1;
                        tp[0][1] = j;
                        tp[0][2] = k;
                        tp[1][0] = i + 1;
                        tp[1][1] = j;
                        tp[1][2] = k + 1;
                        tp[2][0] = i;
                        tp[2][1] = j;
                        tp[2][2] = k + 1;

                    } else if (vp101 && vp001 && vp000) {
                        tp[0][0] = i + 1;
                        tp[0][1] = j;
                        tp[0][2] = k + 1;
                        tp[1][0] = i;
                        tp[1][1] = j;
                        tp[1][2] = k + 1;
                        tp[2][0] = i;
                        tp[2][1] = j;
                        tp[2][2] = k;
                    } else if (vp001 && vp000 && vp100) {
                        tp[0][0] = i;
                        tp[0][1] = j;
                        tp[0][2] = k + 1;
                        tp[1][0] = i;
                        tp[1][1] = j;
                        tp[1][2] = k;
                        tp[2][0] = i + 1;
                        tp[2][1] = j;
                        tp[2][2] = k;
                    }
                    for (ii = 0; ii < 3; ii++) {
                        if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                            vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                            verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                            
                        }
                    }
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]],  vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                    
                }
            }
        }
    }
    console.log(3);
    //face5
    for (i = 0; i < pLength - 1; i++) {
        for (j = 0; j < pWidth - 1; j++) {
            for (k = 0; k < 1; k++) {
                vp000 = vp[pWidth * pHeight * i + pHeight * j + k].isdone;
                vp010 = vp[pWidth * pHeight * i + pHeight * (j + 1) + k].isdone;
                vp100 = vp[pWidth * pHeight * (i + 1) + pHeight * j + k].isdone;
                vp110 = vp[pWidth * pHeight * (i + 1) + pHeight * (j + 1) + k].isdone;

                if (vp000 && vp100 && vp110 && vp010) {
                    tp[0][0] = i;
                    tp[0][1] = j;
                    tp[0][2] = k;
                    tp[1][0] = i + 1;
                    tp[1][1] = j;
                    tp[1][2] = k;
                    tp[2][0] = i + 1;
                    tp[2][1] = j + 1;
                    tp[2][2] = k;
                    tp[3][0] = i;
                    tp[3][1] = j + 1;
                    tp[3][2] = k;

                    for (ii = 0; ii < 4; ii++) {
                        if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                            vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                            verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                            
                        }
                    }
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]],  vertseq[tp[1][0]][tp[1][1]][tp[1][2]]));
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[3][0]][tp[3][1]][tp[3][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                } else if ((vp000 && vp100 && vp110)
                        || (vp100 && vp110 && vp010)
                        || (vp110 && vp010 && vp000)
                        || (vp010 && vp000 && vp100)) {
                    if (vp000 && vp100 && vp110) {
                        tp[0][0] = i;
                        tp[0][1] = j;
                        tp[0][2] = k;
                        tp[1][0] = i + 1;
                        tp[1][1] = j;
                        tp[1][2] = k;
                        tp[2][0] = i + 1;
                        tp[2][1] = j + 1;
                        tp[2][2] = k;

                    } else if (vp100 && vp110 && vp010) {
                        tp[0][0] = i + 1;
                        tp[0][1] = j;
                        tp[0][2] = k;
                        tp[1][0] = i + 1;
                        tp[1][1] = j + 1;
                        tp[1][2] = k;
                        tp[2][0] = i;
                        tp[2][1] = j + 1;
                        tp[2][2] = k;
                    } else if (vp110 && vp010 && vp000) {
                        tp[0][0] = i + 1;
                        tp[0][1] = j + 1;
                        tp[0][2] = k;
                        tp[1][0] = i;
                        tp[1][1] = j + 1;
                        tp[1][2] = k;
                        tp[2][0] = i;
                        tp[2][1] = j;
                        tp[2][2] = k;
                    } else if (vp010 && vp000 && vp100) {
                        tp[0][0] = i;
                        tp[0][1] = j + 1;
                        tp[0][2] = k;
                        tp[1][0] = i;
                        tp[1][1] = j;
                        tp[1][2] = k;
                        tp[2][0] = i + 1;
                        tp[2][1] = j;
                        tp[2][2] = k;
                    }
                    for (ii = 0; ii < 3; ii++) {
                        if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                            vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                            verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                            
                        }
                    }
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]],  vertseq[tp[1][0]][tp[1][1]][tp[1][2]]));
                }
            }
        }
    }
    console.log(5);
	//face2
    for (i = pLength - 1; i < pLength; i++) {
        for (j = 0; j < pWidth - 1; j++) {
            for (k = 0; k < pHeight - 1; k++) {
                vp000 = vp[pWidth * pHeight * i + pHeight * j + k].isdone;
                vp001 = vp[pWidth * pHeight * i + pHeight * j + k + 1].isdone;
                vp010 = vp[pWidth * pHeight * i + pHeight * (j + 1) + k].isdone;
                vp011 = vp[pWidth * pHeight * i + pHeight * (j + 1) + k + 1].isdone;

                if (vp000 && vp010 && vp011 && vp001) {
                    tp[0][0] = i;
                    tp[0][1] = j;
                    tp[0][2] = k;
                    tp[1][0] = i;
                    tp[1][1] = j + 1;
                    tp[1][2] = k;
                    tp[2][0] = i;
                    tp[2][1] = j + 1;
                    tp[2][2] = k + 1;
                    tp[3][0] = i;
                    tp[3][1] = j;
                    tp[3][2] = k + 1;
                    for (ii = 0; ii < 4; ii++) {
                        if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                            vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                            verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                            
                        }
                    }
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]],  vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]],  vertseq[tp[3][0]][tp[3][1]][tp[3][2]]));
                } else if ((vp000 && vp010 && vp011)
                            || (vp010 && vp011 && vp001)
                            || (vp011 && vp001 && vp000)
                            || (vp001 && vp000 && vp010)) {
                    if (vp000 && vp010 && vp011) {
                        tp[0][0] = i;
                        tp[0][1] = j;
                        tp[0][2] = k;
                        tp[1][0] = i;
                        tp[1][1] = j + 1;
                        tp[1][2] = k;
                        tp[2][0] = i;
                        tp[2][1] = j + 1;
                        tp[2][2] = k + 1;
                    } else if (vp010 && vp011 && vp001) {
                        tp[0][0] = i;
                        tp[0][1] = j + 1;
                        tp[0][2] = k;
                        tp[1][0] = i;
                        tp[1][1] = j + 1;
                        tp[1][2] = k + 1;
                        tp[2][0] = i;
                        tp[2][1] = j;
                        tp[2][2] = k + 1;
                    } else if (vp011 && vp001 && vp000) {
                        tp[0][0] = i;
                        tp[0][1] = j + 1;
                        tp[0][2] = k + 1;
                        tp[1][0] = i;
                        tp[1][1] = j;
                        tp[1][2] = k + 1;
                        tp[2][0] = i;
                        tp[2][1] = j;
                        tp[2][2] = k;
                    } else if (vp001 && vp000 && vp010) {
                        tp[0][0] = i;
                        tp[0][1] = j;
                        tp[0][2] = k + 1;
                        tp[1][0] = i;
                        tp[1][1] = j;
                        tp[1][2] = k;
                        tp[2][0] = i;
                        tp[2][1] = j + 1;
                        tp[2][2] = k;
                    }
                    for (ii = 0; ii < 3; ii++) {
                        if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                            vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                            verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                            
                        }
                    }
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]],  vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                }
            }
        }
    }
    console.log(2);
    //face4
    for (i = 0; i < pLength - 1; i++) {
        for (j = pWidth - 1; j < pWidth; j++) {
            for (k = 0; k < pHeight - 1; k++) {
                vp000 = vp[pWidth * pHeight * i + pHeight * j + k].isdone;
                vp001 = vp[pWidth * pHeight * i + pHeight * j + k + 1].isdone;
                vp100 = vp[pWidth * pHeight * (i + 1) + pHeight * j + k].isdone;
                vp101 = vp[pWidth * pHeight * (i + 1) + pHeight * j + k + 1].isdone;

                if (vp000 && vp100 && vp101 && vp001) {
                    tp[0][0] = i;
                    tp[0][1] = j;
                    tp[0][2] = k;
                    tp[1][0] = i + 1;
                    tp[1][1] = j;
                    tp[1][2] = k;
                    tp[2][0] = i + 1;
                    tp[2][1] = j;
                    tp[2][2] = k + 1;
                    tp[3][0] = i;
                    tp[3][1] = j;
                    tp[3][2] = k + 1;
                    for (ii = 0; ii < 4; ii++) {
                        if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                            vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                            verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                            
                        }
                    }
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]]));
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[3][0]][tp[3][1]][tp[3][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                } else if ((vp000 && vp100 && vp101) || (vp100 && vp101 && vp001) || (vp101 && vp001 && vp000) || (vp001 && vp000 && vp100)) {
                    if (vp000 && vp100 && vp101) {
                        tp[0][0] = i;
                        tp[0][1] = j;
                        tp[0][2] = k;
                        tp[1][0] = i + 1;
                        tp[1][1] = j;
                        tp[1][2] = k;
                        tp[2][0] = i + 1;
                        tp[2][1] = j;
                        tp[2][2] = k + 1;
                    } else if (vp100 && vp101 && vp001) {
                        tp[0][0] = i + 1;
                        tp[0][1] = j;
                        tp[0][2] = k;
                        tp[1][0] = i + 1;
                        tp[1][1] = j;
                        tp[1][2] = k + 1;
                        tp[2][0] = i;
                        tp[2][1] = j;
                        tp[2][2] = k + 1;
                    } else if (vp101 && vp001 && vp000) {
                        tp[0][0] = i + 1;
                        tp[0][1] = j;
                        tp[0][2] = k + 1;
                        tp[1][0] = i;
                        tp[1][1] = j;
                        tp[1][2] = k + 1;
                        tp[2][0] = i;
                        tp[2][1] = j;
                        tp[2][2] = k;
                    } else if (vp001 && vp000 && vp100) {
                        tp[0][0] = i;
                        tp[0][1] = j;
                        tp[0][2] = k + 1;
                        tp[1][0] = i;
                        tp[1][1] = j;
                        tp[1][2] = k;
                        tp[2][0] = i + 1;
                        tp[2][1] = j;
                        tp[2][2] = k;
                    }
                    for (ii = 0; ii < 3; ii++) {
                        if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                            vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                            verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                            
                        }
                    }
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]]));
                }

            }
        }
    }
    console.log(4);
    //face6
    for (i = 0; i < pLength - 1; i++) {
        for (j = 0; j < pWidth - 1; j++) {
            for (k = pHeight - 1; k < pHeight; k++) {
                vp000 = vp[pWidth * pHeight * i + pHeight * j + k].isdone;
                vp010 = vp[pWidth * pHeight * i + pHeight * (j + 1) + k].isdone;
                vp100 = vp[pWidth * pHeight * (i + 1) + pHeight * j + k].isdone;
                vp110 = vp[pWidth * pHeight * (i + 1) + pHeight * (j + 1) + k].isdone;

                if (vp000 && vp100 && vp110 && vp010) {
                    tp[0][0] = i;
                    tp[0][1] = j;
                    tp[0][2] = k;
                    tp[1][0] = i + 1;
                    tp[1][1] = j;
                    tp[1][2] = k;
                    tp[2][0] = i + 1;
                    tp[2][1] = j + 1;
                    tp[2][2] = k;
                    tp[3][0] = i;
                    tp[3][1] = j + 1;
                    tp[3][2] = k;
                    for (ii = 0; ii < 4; ii++) {
                        if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                            vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                            verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                            
                        }
                    }
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]], vertseq[tp[3][0]][tp[3][1]][tp[3][2]]));
                } else if ((vp000 && vp100 && vp110)
                               || (vp100 && vp110 && vp010)
                               || (vp110 && vp010 && vp000)
                               || (vp010 && vp000 && vp100)) {
                    if (vp000 && vp100 && vp110) {
                        tp[0][0] = i;
                        tp[0][1] = j;
                        tp[0][2] = k;
                        tp[1][0] = i + 1;
                        tp[1][1] = j;
                        tp[1][2] = k;
                        tp[2][0] = i + 1;
                        tp[2][1] = j + 1;
                        tp[2][2] = k;
                    } else if (vp100 && vp110 && vp010) {
                        tp[0][0] = i + 1;
                        tp[0][1] = j;
                        tp[0][2] = k;
                        tp[1][0] = i + 1;
                        tp[1][1] = j + 1;
                        tp[1][2] = k;
                        tp[2][0] = i;
                        tp[2][1] = j + 1;
                        tp[2][2] = k;
                    } else if (vp110 && vp010 && vp000) {
                        tp[0][0] = i + 1;
                        tp[0][1] = j + 1;
                        tp[0][2] = k;
                        tp[1][0] = i;
                        tp[1][1] = j + 1;
                        tp[1][2] = k;
                        tp[2][0] = i;
                        tp[2][1] = j;
                        tp[2][2] = k;
                    } else if (vp010 && vp000 && vp100) {
                        tp[0][0] = i;
                        tp[0][1] = j + 1;
                        tp[0][2] = k;
                        tp[1][0] = i;
                        tp[1][1] = j;
                        tp[1][2] = k;
                        tp[2][0] = i + 1;
                        tp[2][1] = j;
                        tp[2][2] = k;
                    }
                    for (ii = 0; ii < 3; ii++) {
                        if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                            vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                            verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                            
                        }
                    }
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                    
                }
            }
        }
    }
    console.log(6);
    for (i = 0; i < pLength - 1; i++) {
        console.log(i);
        for (j = 0; j < pWidth - 1; j++) {
            for (k = 0; k < pHeight - 1; k++) {
                vp000 = vp[pWidth * pHeight * i + pHeight * j + k].isdone;
                vp001 = vp[pWidth * pHeight * i + pHeight * j + k + 1].isdone;
                vp010 = vp[pWidth * pHeight * i + pHeight * (j + 1) + k].isdone;
                vp011 = vp[pWidth * pHeight * i + pHeight * (j + 1) + k + 1].isdone;
                vp100 = vp[pWidth * pHeight * (i + 1) + pHeight * j + k].isdone;
                vp101 = vp[pWidth * pHeight * (i + 1) + pHeight * j + k + 1].isdone;
                vp110 = vp[pWidth * pHeight * (i + 1) + pHeight * (j + 1) + k].isdone;
                vp111 = vp[pWidth * pHeight * (i + 1) + pHeight * (j + 1) + k + 1].isdone;


                sumtype = 0;
                for (ii = 0; ii < 2; ii++) {
                    for (jj = 0; jj < 2; jj++) {
                        for (kk = 0; kk < 2; kk++) {
                            if (vp[pWidth * pHeight * (i + ii) + pHeight * (j + jj) + k + kk].isdone) {
                                sumtype++;
                            }
                        }
                    }
                }

                if (sumtype === 3) {
                    if ((vp000 && vp100 && vp110)
                            || (vp000 && vp010 && vp110)
                            || (vp010 && vp100 && vp110)
                            || (vp000 && vp010 && vp100)
                            || (vp001 && vp101 && vp111)
                            || (vp001 && vp011 && vp111)
                            || (vp011 && vp101 && vp111)
                            || (vp001 && vp011 && vp101)
                            || (vp000 && vp100 && vp101)
                            || (vp100 && vp101 && vp001)
                            || (vp000 && vp101 && vp001)
                            || (vp000 && vp100 && vp001)
                            || (vp110 && vp100 && vp111)
                            || (vp110 && vp101 && vp111)
                            || (vp100 && vp101 && vp111)
                            || (vp110 && vp100 && vp101)
                            || (vp110 && vp010 && vp011)
                            || (vp010 && vp011 && vp111)
                            || (vp110 && vp011 && vp111)
                            || (vp110 && vp010 && vp111)
                            || (vp000 && vp010 && vp001)
                            || (vp000 && vp001 && vp011)
                            || (vp001 && vp010 && vp011)
                            || (vp000 && vp010 && vp011)) {
                        if (vp000 && vp100 && vp110) {
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                        } else if (vp000 && vp010 && vp110) { // 11
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;
                        } else if (vp010 && vp100 && vp110) { // 12
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                        } else if (vp000 && vp010 && vp100) { // 13
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;
                        } else if (vp001 && vp101 && vp111) { // 14
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                        } else if (vp001 && vp011 && vp111) { // 21
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                        } else if (vp011 && vp101 && vp111) { // 22
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                        } else if (vp001 && vp011 && vp101) { // 23
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                        } else if (vp000 && vp100 && vp101) { // 24
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;
                        } else if (vp100 && vp101 && vp001) { // 31
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;
                        } else if (vp000 && vp101 && vp001) { // 32
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                        } else if (vp000 && vp100 && vp001) { // 33
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                        } else if (vp110 && vp100 && vp111) { // 34
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;
                        } else if (vp110 && vp101 && vp111) { // 41
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                        } else if (vp100 && vp101 && vp111) { // 42
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                        } else if (vp110 && vp100	&& vp101) { // 43
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                        } else if (vp110 && vp010 && vp011) { // 44
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                        } else if (vp010 && vp011 && vp111) { // 51
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                        } else if (vp110 && vp011 && vp111) { // 52
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                        } else if (vp110 && vp010 && vp111) { // 53
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                        } else if (vp000 && vp010	&& vp001) { // 54
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                        } else if (vp000 && vp001 && vp011) { // 61
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;
                        } else if (vp001 && vp010 && vp011) { // 62
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                        } else if (vp000 && vp010 && vp011) { // 62
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                        }//64
                        for (ii = 0; ii < 3; ii++) {
                            if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                                vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                                verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                                
                            }
                        }
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                        
                    }//no5 24
                //total3
                } else if (sumtype === 4) { // CHECK
                    if ((vp000 && vp100 && vp110 && vp010)
                            || (vp001 && vp101 && vp111 && vp011)
                            || (vp000 && vp100 && vp101 && vp001)
                            || (vp110 && vp100 && vp101 && vp111)
                            || (vp110 && vp010 && vp011 && vp111)
                            || (vp000 && vp010 && vp001 && vp011)) {
                        if (vp000 && vp100 && vp110 && vp010) {
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;
                        } else if (vp001 && vp101 && vp111 && vp011) {
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                        } else if (vp000 && vp100 && vp101 && vp001) {
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                        } else if (vp110 && vp100 && vp101 && vp111) {
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                        } else if (vp110 && vp010 && vp011 && vp111) {
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                        } else if (vp000 && vp010 && vp001 && vp011) {
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;
                        }
                        for (ii = 0; ii < 4; ii++) {
                            if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                                vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                                verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                                
                            }
                        }
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                        
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]], vertseq[tp[3][0]][tp[3][1]][tp[3][2]]));
                        
                    //no.8 6
                    } else if ((vp000 && vp100 && vp110  && vp011)//11
                                 || (vp000 && vp010 && vp110 && vp101)//12
                                 || (vp010 && vp100 && vp110 && vp001)//13
                                 || (vp000 && vp010 && vp100 && vp111)//14
                                 || (vp001 && vp101 && vp111 && vp010)//21
                                 || (vp001 && vp011 && vp111 && vp100)//22
                                 || (vp011 && vp101 && vp111 && vp000)//23
                                 || (vp001 && vp011 && vp101 && vp110)//24
                                 || (vp000 && vp100 && vp101 && vp011)//31
                                 || (vp100 && vp101 && vp001 && vp010)//32
                                 || (vp000 && vp101 && vp001 && vp110)//33
                                 || (vp000 && vp100 && vp001 && vp111)//34
                                 || (vp110 && vp100 && vp111 && vp001)//41
                                 || (vp110 && vp101 && vp111 && vp000)//42
                                 || (vp100 && vp101 && vp111 && vp010)//43
                                 || (vp110 && vp100 && vp101 && vp011)//44
                                 || (vp110 && vp010 && vp011 && vp101)//51
                                 || (vp010 && vp011 && vp111 && vp100)//52
                                 || (vp110 && vp011 && vp111 && vp000)//53
                                 || (vp110 && vp010 && vp111 && vp001)//54
                                 || (vp000 && vp010 && vp001 && vp111)//61
                                 || (vp000 && vp001 && vp011 && vp110)//62
                                 || (vp001 && vp010 && vp011 && vp100)//63
                                 || (vp000 && vp010 && vp011 && vp101)) {
                        if (vp000 && vp100 && vp110  && vp011) {
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                        } else if (vp000 && vp010 && vp110 && vp101) { // 11
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;
                        } else if (vp010 && vp100 && vp110 && vp001) { //12
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                        } else if (vp000 && vp010 && vp100 && vp111) { //13
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;
                        } else if (vp001 && vp101 && vp111 && vp010) { //14
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                        } else if (vp001 && vp011 && vp111 && vp100) { //21
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                        } else if (vp011 && vp101 && vp111 && vp000) { //22
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                        } else if (vp001 && vp011 && vp101 && vp110) { //23
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                        } else if (vp000 && vp100 && vp101 && vp011) { //24
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;
                        } else if (vp100 && vp101 && vp001 && vp010) { //31
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;
                        } else if (vp000 && vp101 && vp001 && vp110) { //32
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                        } else if (vp000 && vp100 && vp001 && vp111) { //33
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                        } else if (vp110 && vp100 && vp111 && vp001) { //34
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;
                        } else if (vp110 && vp101 && vp111 && vp000) { //41
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                        } else if (vp100 && vp101 && vp111 && vp010) { //42
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                        } else if (vp110 && vp100	&& vp101 && vp011) { //43
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                        } else if (vp110 && vp010 && vp011 && vp101) { //44
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                        } else if (vp010 && vp011 && vp111 && vp100) { //51
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                        } else if (vp110 && vp011 && vp111 && vp000) { //52
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                        } else if (vp110 && vp010 && vp111 && vp001) { //53
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                        } else if (vp000 && vp010	&& vp001 && vp111) { //54
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                        } else if (vp000 && vp001 && vp011 && vp110) { //61
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;
                        } else if (vp001 && vp010 && vp011 && vp100) { //62
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                        } else if (vp000 && vp010	&& vp011 && vp101) { //63
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                        }
                        for (ii = 0; ii < 3; ii++) { //64
                            if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                                vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                                verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                                
                            }
                        }
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                        
                    } else if ((vp000 && vp011 && vp110 && vp010)//no12 24
                                || (vp000 && vp100 && vp110 && vp101)
                                || (vp000 && vp001 && vp100 && vp010)
                                || (vp010 && vp100 && vp110 && vp111)
                                || (vp001 && vp011 && vp111 && vp010)
                                || (vp001 && vp100 && vp111 && vp101)
                                || (vp000 && vp001 && vp101 && vp011)
                                || (vp011 && vp101 && vp110 && vp111)) {
                        if (vp010 && vp011 && vp000 && vp110) {
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                        } else if (vp100 && vp101 && vp110 && vp000) { //1
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;
                        } else if (vp000 && vp001 && vp100 && vp010) { //2
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                        } else if (vp110 && vp111 && vp010 && vp100) { //3
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;
                        } else if (vp010 && vp011 && vp111 && vp001) { //4
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                        } else if (vp100 && vp101 && vp111 && vp001) { //5
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                        } else if (vp000 && vp001 && vp101 && vp011) { //6
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                        } else if (vp011 && vp101 && vp110 && vp111) { //7
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                        }//8
                        for (ii = 0; ii < 3; ii++) {
                            if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                                vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                                verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                                
                            }
                        }
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                        
                    // no.9 8
                    } else if ((vp000 && vp100 && vp110 && vp001)
                                 || (vp010 && vp100 && vp110 && vp101)
                                 || (vp010 && vp000 && vp110 && vp111)
                                 || (vp010 && vp000 && vp100 && vp011)
                                 || (vp011 && vp001 && vp101 && vp100)
                                 || (vp111 && vp001 && vp101 && vp110)
                                 || (vp111 && vp011 && vp101 && vp010)
                                 || (vp111 && vp011 && vp001 && vp000)
                                 || (vp110 && vp011 && vp001 && vp010)
                                 || (vp101 && vp000 && vp001 && vp010)
                                 || (vp101 && vp000 && vp111 && vp100)
                                 || (vp011 && vp110 && vp111 && vp100)) {
                        if (vp000 && vp100 && vp110 && vp001) {
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;
                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;
                        } else if (vp010 && vp100 && vp110 && vp101) { //1
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;
                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;
                        } else if (vp010 && vp000 && vp110 && vp111) { //2
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                        } else if (vp010 && vp000 && vp100 && vp011) { //3
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                        } else if (vp011 && vp001 && vp101 && vp100) { //4
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                        } else if (vp111 && vp001 && vp101 && vp110) { //5
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;
                        } else if (vp111 && vp011 && vp101 && vp010) { //6
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;
                        } else if (vp111 && vp011 && vp001 && vp000) { //7
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                        } else if (vp110 && vp011 && vp001 && vp010) { //8
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;
                        } else if (vp101 && vp000 && vp001 && vp010) { //9
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;
                        } else if (vp101 && vp000 && vp111 && vp100) { //10
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                        } else if (vp011 && vp110 && vp111 && vp100) { //11
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                        } //12
                        for (ii = 0; ii < 4; ii++) {
                            if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                                vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                                verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                                
                            }
                        }
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                        
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]], vertseq[tp[3][0]][tp[3][1]][tp[3][2]]));
                        
                    //no.11 12
                    } else if ((vp000 && vp100 && vp010 && vp101)
                                 || (vp000 && vp100 && vp110 && vp111)
                                 || (vp010 && vp100 && vp110 && vp011)
                                 || (vp010 && vp000 && vp110 && vp001)
                                 || (vp111 && vp001 && vp101 && vp000)
                                 || (vp111 && vp011 && vp101 && vp100)
                                 || (vp111 && vp011 && vp001 && vp110)
                                 || (vp101 && vp011 && vp001 && vp010)
                                 || (vp111 && vp011 && vp000 && vp010)
                                 || (vp100 && vp000 && vp001 && vp011)
                                 || (vp101 && vp001 && vp110 && vp100)
                                 || (vp010 && vp110 && vp111 && vp101)) {
                        if (vp000 && vp100 && vp010 && vp101) {
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;
                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;
                        } else if (vp000 && vp100 && vp110 && vp111) { //1
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;
                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                        } else if (vp010 && vp100 && vp110 && vp011) { //2
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                        } else if (vp010 && vp000 && vp110 && vp001) { //3
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;
                        } else if (vp111 && vp001 && vp101 && vp000) { //4
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                        } else if (vp111 && vp011 && vp101 && vp100) { //5
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                        } else if (vp111 && vp011 && vp001 && vp110) { //6
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;
                        } else if (vp101 && vp011 && vp001 && vp010) { //7
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;
                        } else if (vp111 && vp011 && vp000 && vp010) { //8
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k;
                        } else if (vp100 && vp000 && vp001 && vp011) { //9
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k;
                        } else if (vp101 && vp001 && vp110 && vp100) { //10
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;
                        } else if (vp010 && vp110 && vp111 && vp101) { //11
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;
                        }//12
                        for (ii = 0; ii < 4; ii++) {
                            if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                                vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                                verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                                
                            }
                        }
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                        
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]], vertseq[tp[3][0]][tp[3][1]][tp[3][2]]));
                        
                    }//no.14 12
                    //total4
                } else if (sumtype === 5) {
                    if ((!vp100 && !vp001 && !vp111)
                            || (!vp010 && !vp001 && !vp111)
                            || (!vp110 && !vp101 && !vp011)
                            || (!vp000 && !vp101 && !vp011)
                            || (!vp101 && !vp000 && !vp110)
                            || (!vp011 && !vp000 && !vp110)
                            || (!vp111 && !vp100 && !vp010)
                            || (!vp001 && !vp100 && !vp010)) {
                        if (!vp100 && !vp001 && !vp111) {
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                        } else if (!vp010 && !vp001 && !vp111) { //1
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;

                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;

                        } else if (!vp110 && !vp101 && !vp011) { //2
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                        } else if (!vp000 && !vp101 && !vp011) { //3
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;

                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;

                        } else if (!vp101 && !vp000 && !vp110) { //4
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;

                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                        } else if (!vp011 && !vp000 && !vp110) { //5
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;

                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                        } else if (!vp111 && !vp100 && !vp010) { //6
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;

                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                        } else if (!vp001 && !vp100 && !vp010) { //7
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;

                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                        }//8
                        for (ii = 0; ii < 3; ii++) {
                            if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                                vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                                verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                                
                            }
                        }
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                        
                    //no.7 8
                    } else if ((!vp000 && !vp100 && !vp110)
                             || (!vp000 && !vp010 && !vp110)
                             || (!vp010 && !vp100 && !vp110)
                             || (!vp000 && !vp010 && !vp100)
                             || (!vp001 && !vp101 && !vp111)
                             || (!vp001 && !vp011 && !vp111)
                             || (!vp011 && !vp101 && !vp111)
                             || (!vp001 && !vp011 && !vp101)
                             || (!vp000 && !vp100 && !vp101)
                             || (!vp100 && !vp101 && !vp001)
                             || (!vp000 && !vp101 && !vp001)
                             || (!vp000 && !vp100 && !vp001)
                             || (!vp110 && !vp100 && !vp111)
                             || (!vp110 && !vp101 && !vp111)
                             || (!vp100 && !vp101 && !vp111)
                             || (!vp110 && !vp100 && !vp101)
                             || (!vp110 && !vp010 && !vp011)
                             || (!vp010 && !vp011 && !vp111)
                             || (!vp110 && !vp011 && !vp111)
                             || (!vp110 && !vp010 && !vp111)
                             || (!vp000 && !vp010 && !vp001)
                             || (!vp000 && !vp001 && !vp011)
                             || (!vp001 && !vp010 && !vp011)
                             || (!vp000 && !vp010 && !vp011)) {
                        if (!vp000 && !vp100 && !vp110) {
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                        } else if (!vp000 && !vp010 && !vp110) { //11
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k;

                        } else if (!vp010 && !vp100 && !vp110) { //12
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k;

                        } else if (!vp000 && !vp010 && !vp100) { //13
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                        } else if (!vp001 && !vp101 && !vp111) { //14
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                        } else if (!vp001 && !vp011 && !vp111) { //21
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;

                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;

                        } else if (!vp011 && !vp101 && !vp111) { //22
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;

                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;

                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;

                        } else if (!vp001 && !vp011 && !vp101) { //23
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                        } else if (!vp000 && !vp100 && !vp101) { //24
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;

                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;

                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;

                        } else if (!vp100 && !vp101 && !vp001) { //31
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;

                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k;

                        } else if (!vp000 && !vp101 && !vp001) { //32
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k;

                        } else if (!vp000 && !vp100 && !vp001) { //33
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;

                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;

                        } else if (!vp110 && !vp100 && !vp111) { //34
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;

                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;

                        } else if (!vp110 && !vp101 && !vp111) { //41
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k;

                        } else if (!vp100 && !vp101 && !vp111) { //42
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                        } else if (!vp110 && !vp100 && !vp101) { //43
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                        } else if (!vp110 && !vp010 && !vp011) { //44
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                        } else if (!vp010 && !vp011 && !vp111) { //51
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                        } else if (!vp110 && !vp011 && !vp111) { //52
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                        } else if (!vp110 && !vp010 && !vp111) { //53
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                        } else if (!vp000 && !vp010 && !vp001) { //54
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                        } else if (!vp000 && !vp001 && !vp011) { //61
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                        } else if (!vp001 && !vp010 && !vp011) { //62
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k;

                        } else if (!vp000 && !vp010 && !vp011) { //63
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;

                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;

                        }
                        for (ii = 0; ii < 4; ii++) { //64
                            if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                                vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                                verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                                
                            }
                        }
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                        
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]], vertseq[tp[3][0]][tp[3][1]][tp[3][2]]));
                        
                    } else if ((!vp000 && !vp100 && !vp111)//1 //no5 24
                                 || (!vp010 && !vp110 && !vp001)//2
                                 || (!vp011 && !vp111 && !vp100)//3
                                 || (!vp001 && !vp101 && !vp110)//4
                                 || (!vp000 && !vp010 && !vp111)//5
                                 || (!vp101 && !vp111 && !vp010)//6
                                 || (!vp100 && !vp110 && !vp011)//7
                                 || (!vp001 && !vp011 && !vp110)//8
                                 || (!vp000 && !vp001 && !vp111)//9
                                 || (!vp110 && !vp111 && !vp000)//10
                                 || (!vp100 && !vp101 && !vp011)//11
                                 || (!vp010 && !vp011 && !vp101)) {
                        if (!vp000 && !vp100 && !vp111) {
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;
                            tp[4][0] = i;
                            tp[4][1] = j + 1;
                            tp[4][2] = k + 1;

                        } else if (!vp010 && !vp110 && !vp001) { //1
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;
                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                            tp[4][0] = i + 1;
                            tp[4][1] = j;
                            tp[4][2] = k + 1;
                        } else if (!vp011 && !vp111 && !vp100) { //2
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;
                            tp[4][0] = i;
                            tp[4][1] = j;
                            tp[4][2] = k;
                        } else if (!vp001 && !vp101 && !vp110) { //3
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;
                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                            tp[4][0] = i;
                            tp[4][1] = j + 1;
                            tp[4][2] = k;
                        } else if (!vp000 && !vp010 && !vp111) { //4
                                 //tp[0][0] = i;tp[0][1] = j;tp[0][2] = k + 1;
                                 //tp[1][0] = i + 1;tp[1][1] = j;tp[1][2] = k;
                                 //tp[2][0] = i + 1;tp[2][1] = j + 1;tp[2][2] = k;
                                 //tp[3][0] = i;tp[3][1] = j + 1;tp[3][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[4][0] = i + 1;
                            tp[4][1] = j;
                            tp[4][2] = k + 1;
                        } else if (!vp101 && !vp111 && !vp010) { //5
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                            tp[4][0] = i;
                            tp[4][1] = j;
                            tp[4][2] = k;
                        } else if (!vp100 && !vp110 && !vp011) { //6
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                            tp[4][0] = i;
                            tp[4][1] = j;
                            tp[4][2] = k + 1;
                        } else if (!vp001 && !vp011 && !vp110) { //7
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;
                            tp[4][0] = i + 1;
                            tp[4][1] = j;
                            tp[4][2] = k;
                        } else if (!vp000 && !vp001 && !vp111) { // 8
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                            tp[4][0] = i + 1;
                            tp[4][1] = j + 1;
                            tp[4][2] = k;
                        } else if (!vp110 && !vp111 && !vp000) { // 9
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;
                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;
                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;
                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;
                            tp[4][0] = i;
                            tp[4][1] = j;
                            tp[4][2] = k + 1;
                        } else if (!vp100 && !vp101 && !vp011) { //10
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;
                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;
                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;
                            tp[4][0] = i;
                            tp[4][1] = j + 1;
                            tp[4][2] = k;
                        } else if (!vp010 && !vp011 && !vp101) { //11
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;
                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;
                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;
                            tp[4][0] = i + 1;
                            tp[4][1] = j;
                            tp[4][2] = k;
                        }//12
                        for (ii = 0; ii < 5; ii++) {
                            if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                                vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                                verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                                
                            }
                        }
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                        
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]], vertseq[tp[3][0]][tp[3][1]][tp[3][2]]));
                        
                        faces.push(new THREE.Face3(vertseq[tp[2][0]][tp[2][1]][tp[2][2]], vertseq[tp[4][0]][tp[4][1]][tp[4][2]], vertseq[tp[3][0]][tp[3][1]][tp[3][2]]));
                        
                    //no.6 12-1
                    } else if ((!vp000 && !vp100 && !vp011)//1
                             || (!vp010 && !vp110 && !vp101)//2
                             || (!vp011 && !vp111 && !vp000)//3
                             || (!vp001 && !vp101 && !vp010)//4
                             || (!vp000 && !vp010 && !vp101)//5
                             || (!vp101 && !vp111 && !vp000)//6
                             || (!vp100 && !vp110 && !vp001)//7
                             || (!vp001 && !vp011 && !vp100)//8
                             || (!vp000 && !vp001 && !vp110)//9
                             || (!vp110 && !vp111 && !vp001)//10
                             || (!vp100 && !vp101 && !vp010)//11
                             || (!vp010 && !vp011 && !vp100)) {
                        if (!vp000 && !vp100 && !vp011) {
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;
                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;
                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;
                            tp[4][0] = i + 1;
                            tp[4][1] = j + 1;
                            tp[4][2] = k + 1;
                        } else if (!vp010 && !vp110 && !vp101) { //1
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                            tp[4][0] = i;
                            tp[4][1] = j;
                            tp[4][2] = k + 1;

                        } else if (!vp011 && !vp111 && !vp000) { //2
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;

                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;

                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                            tp[4][0] = i + 1;
                            tp[4][1] = j;
                            tp[4][2] = k;

                        } else if (!vp001 && !vp101 && !vp010) { //3
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                            tp[4][0] = i + 1;
                            tp[4][1] = j + 1;
                            tp[4][2] = k;

                        } else if (!vp000 && !vp010 && !vp101) { //4
                                //tp[0][0] = i;tp[0][1] = j;tp[0][2] = k + 1;
                                //tp[1][0] = i + 1;tp[1][1] = j;tp[1][2] = k;
                                //tp[2][0] = i + 1;tp[2][1] = j + 1;tp[2][2] = k;
                                //tp[3][0] = i;tp[3][1] = j + 1;tp[3][2] = k + 1;
                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;

                            tp[4][0] = i + 1;
                            tp[4][1] = j + 1;
                            tp[4][2] = k + 1;

                        } else if (!vp101 && !vp111 && !vp000) { //5
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                            tp[4][0] = i;
                            tp[4][1] = j + 1;
                            tp[4][2] = k;

                        } else if (!vp100 && !vp110 && !vp001) { //6
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                            tp[4][0] = i;
                            tp[4][1] = j + 1;
                            tp[4][2] = k + 1;

                        } else if (!vp001 && !vp011 && !vp100) { //7
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                            tp[4][0] = i + 1;
                            tp[4][1] = j + 1;
                            tp[4][2] = k;

                        } else if (!vp000 && !vp001 && !vp110) { //8
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                            tp[4][0] = i + 1;
                            tp[4][1] = j + 1;
                            tp[4][2] = k + 1;

                        } else if (!vp110 && !vp111 && !vp001) { //9
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                            tp[4][0] = i;
                            tp[4][1] = j;
                            tp[4][2] = k;

                        } else if (!vp100 && !vp101 && !vp010) { //10
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;

                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;

                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;

                            tp[4][0] = i;
                            tp[4][1] = j + 1;
                            tp[4][2] = k + 1;

                        } else if (!vp010 && !vp011 && !vp100) { //11
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;

                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;

                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                            tp[4][0] = i + 1;
                            tp[4][1] = j;
                            tp[4][2] = k + 1;

                        }//12
                        for (ii = 0; ii < 5; ii++) {
                            if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                                vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                                verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                                
                            }
                        }
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                        
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]], vertseq[tp[3][0]][tp[3][1]][tp[3][2]]));
                        
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[4][0]][tp[4][1]][tp[4][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]]));
                        
                    }//no.6 12-2
                //total5
                } else if (sumtype === 6) {
                    if ((!vp000 && !vp100)
                            || (!vp010 && !vp110)
                            || (!vp011 && !vp111)
                            || (!vp001 && !vp101)
                            || (!vp000 && !vp010)
                            || (!vp101 && !vp111)
                            || (!vp100 && !vp110)
                            || (!vp001 && !vp011)
                            || (!vp000 && !vp001)
                            || (!vp110 && !vp111)
                            || (!vp100 && !vp101)
                            || (!vp010 && !vp011)) {
                        if (!vp000 && !vp100) {
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;

                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;

                        } else if (!vp010 && !vp110) { //1
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k;

                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                        } else if (!vp011 && !vp111) { //2
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                        } else if (!vp001 && !vp101) { //3
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                        } else if (!vp000 && !vp010) { //4
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k;

                        } else if (!vp101 && !vp111) { //5
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;

                        } else if (!vp100 && !vp110) { //6
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                        } else if (!vp001 && !vp011) { //7
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;

                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;

                        } else if (!vp000 && !vp001) { //8
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k;

                        } else if (!vp110 && !vp111) { //9
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;

                        } else if (!vp100 && !vp101) { //10
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                        } else if (!vp010 && !vp011) { //11
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;

                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;

                        }//12
                        for (ii = 0; ii < 4; ii++) {
                            if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                                vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                                verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                                
                            }
                        }
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                        
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]], vertseq[tp[3][0]][tp[3][1]][tp[3][2]]));
                        
                    //no.2 12	
                    } else if ((!vp000 && !vp111)
                                 || (!vp100 && !vp011)
                                 || (!vp010 && !vp101)
                                 || (!vp110 && !vp001)) {
                        if (!vp000 && !vp111) {
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;

                            tp[4][0] = i;
                            tp[4][1] = j + 1;
                            tp[4][2] = k;

                            tp[5][0] = i + 1;
                            tp[5][1] = j;
                            tp[5][2] = k;

                        } else if (!vp100 && !vp011) { //1
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;

                            tp[4][0] = i;
                            tp[4][1] = j;
                            tp[4][2] = k;

                            tp[5][0] = i + 1;
                            tp[5][1] = j + 1;
                            tp[5][2] = k;

                        } else if (!vp010 && !vp101) { //2
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                            tp[4][0] = i + 1;
                            tp[4][1] = j + 1;
                            tp[4][2] = k;

                            tp[5][0] = i;
                            tp[5][1] = j;
                            tp[5][2] = k;

                        } else if (!vp110 && !vp001) { //3
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                            tp[4][0] = i + 1;
                            tp[4][1] = j;
                            tp[4][2] = k;

                            tp[5][0] = i;
                            tp[5][1] = j + 1;
                            tp[5][2] = k;

                        } // 4
                        for (ii = 0; ii < 6; ii++) {
                            if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                                vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                                verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                                
                            }
                        }
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                        
                        faces.push(new THREE.Face3(vertseq[tp[3][0]][tp[3][1]][tp[3][2]], vertseq[tp[4][0]][tp[4][1]][tp[4][2]], vertseq[tp[5][0]][tp[5][1]][tp[5][2]]));
                        
                    //no.4 4
                    } else if ((!vp000 && !vp101)
                                 || (!vp100 && !vp001)
                                 || (!vp100 && !vp111)
                                 || (!vp110 && !vp101)
                                 || (!vp110 && !vp011)
                                 || (!vp010 && !vp111)
                                 || (!vp010 && !vp001)
                                 || (!vp000 && !vp011)
                                 || (!vp001 && !vp111)
                                 || (!vp101 && !vp011)
                                 || (!vp000 && !vp110)
                                 || (!vp100 && !vp010)) {
                        if (!vp000 && !vp101) {
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                        } else if (!vp100 && !vp001) { //1
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                        } else if (!vp100 && !vp111) { //2
                            tp[0][0] = i + 1;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k + 1;

                        } else if (!vp110 && !vp101) { //3
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                        } else if (!vp110 && !vp011) { //4
                            tp[0][0] = i + 1;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;

                        } else if (!vp010 && !vp111) { //5
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k + 1;

                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k;

                        } else if (!vp010 && !vp001) { //6
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[2][0] = i;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;

                        } else if (!vp000 && !vp011) { //7
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[2][0] = i;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k;

                        } else if (!vp001 && !vp111) { //8
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k + 1;

                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k + 1;

                            tp[1][0] = i;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[3][0] = i + 1;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                        } else if (!vp101 && !vp011) { //9
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k + 1;

                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k + 1;

                            tp[1][0] = i + 1;
                            tp[1][1] = j;
                            tp[1][2] = k;

                            tp[3][0] = i;
                            tp[3][1] = j + 1;
                            tp[3][2] = k;

                        } else if (!vp000 && !vp110) { //10
                            tp[0][0] = i;
                            tp[0][1] = j + 1;
                            tp[0][2] = k;

                            tp[2][0] = i + 1;
                            tp[2][1] = j;
                            tp[2][2] = k;

                            tp[1][0] = i + 1;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                            tp[3][0] = i;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;

                        } else if (!vp100 && !vp010) { //11
                            tp[0][0] = i;
                            tp[0][1] = j;
                            tp[0][2] = k;

                            tp[2][0] = i + 1;
                            tp[2][1] = j + 1;
                            tp[2][2] = k;

                            tp[1][0] = i;
                            tp[1][1] = j + 1;
                            tp[1][2] = k + 1;

                            tp[3][0] = i + 1;
                            tp[3][1] = j;
                            tp[3][2] = k + 1;

                        }//12
                        for (ii = 0; ii < 4; ii++) {
                            if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                                vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                                verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                                
                            }
                        }
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                        
                        faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]], vertseq[tp[3][0]][tp[3][1]][tp[3][2]]));
                        
                    }//no.3 12
                //total6 
                } else if (sumtype === 7) {
                    if (!vp000) {
                        tp[0][0] = i;
                        tp[0][1] = j + 1;
                        tp[0][2] = k;

                        tp[1][0] = i + 1;
                        tp[1][1] = j;
                        tp[1][2] = k;

                        tp[2][0] = i;
                        tp[2][1] = j;
                        tp[2][2] = k + 1;

                    } else if (!vp100) { //1
                        tp[0][0] = i;
                        tp[0][1] = j;
                        tp[0][2] = k;
                        tp[1][0] = i + 1;
                        tp[1][1] = j + 1;
                        tp[1][2] = k;
                        tp[2][0] = i + 1;
                        tp[2][1] = j;
                        tp[2][2] = k + 1;
                    } else if (!vp110) { //2
                        tp[0][0] = i + 1;
                        tp[0][1] = j;
                        tp[0][2] = k;

                        tp[1][0] = i;
                        tp[1][1] = j + 1;
                        tp[1][2] = k;

                        tp[2][0] = i + 1;
                        tp[2][1] = j + 1;
                        tp[2][2] = k + 1;

                    } else if (!vp010) {       //3
                        tp[0][0] = i + 1;
                        tp[0][1] = j + 1;
                        tp[0][2] = k;

                        tp[1][0] = i;
                        tp[1][1] = j;
                        tp[1][2] = k;

                        tp[2][0] = i;
                        tp[2][1] = j + 1;
                        tp[2][2] = k + 1;

                    } else if (!vp001) { //4
                        tp[0][0] = i + 1;
                        tp[0][1] = j;
                        tp[0][2] = k + 1;

                        tp[1][0] = i;
                        tp[1][1] = j + 1;
                        tp[1][2] = k + 1;

                        tp[2][0] = i;
                        tp[2][1] = j;
                        tp[2][2] = k;

                    } else if (!vp101) { //5
                        tp[0][0] = i + 1;
                        tp[0][1] = j + 1;
                        tp[0][2] = k + 1;

                        tp[1][0] = i;
                        tp[1][1] = j;
                        tp[1][2] = k + 1;

                        tp[2][0] = i + 1;
                        tp[2][1] = j;
                        tp[2][2] = k;

                    } else if (!vp111) { //6
                        tp[0][0] = i;
                        tp[0][1] = j + 1;
                        tp[0][2] = k + 1;

                        tp[1][0] = i + 1;
                        tp[1][1] = j;
                        tp[1][2] = k + 1;

                        tp[2][0] = i + 1;
                        tp[2][1] = j + 1;
                        tp[2][2] = k;

                    } else if (!vp011) { //7
                        tp[0][0] = i;
                        tp[0][1] = j;
                        tp[0][2] = k + 1;

                        tp[1][0] = i + 1;
                        tp[1][1] = j + 1;
                        tp[1][2] = k + 1;

                        tp[2][0] = i;
                        tp[2][1] = j + 1;
                        tp[2][2] = k;

                    }//8
                    for (ii = 0; ii < 3; ii++) {
                        if (vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] === -1) {
                            vertseq[tp[ii][0]][tp[ii][1]][tp[ii][2]] = verts.length;
                            verts.push(new THREE.Vector3(tp[ii][0], tp[ii][1], tp[ii][2]));
                            
                        }
                    }
                    faces.push(new THREE.Face3(vertseq[tp[0][0]][tp[0][1]][tp[0][2]], vertseq[tp[1][0]][tp[1][1]][tp[1][2]], vertseq[tp[2][0]][tp[2][1]][tp[2][2]]));
                    
                }//total7
            }//every ijk
        }//j
    }//i
    this.faces = faces;
    this.verts = verts;
    for (i = 0; i < verts.length; i++) {
        this.verts[i].atomid = this.vp[verts[i].x * pWidth * pHeight + pHeight * verts[i].y + verts[i].z].atomid;
    }
};
