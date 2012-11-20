/*global jQuery, $, console, THREE */

/*
 GLmol - Molecular Viewer on WebGL/Javascript (0.48/G)
  (C) Copyright 2011-2012, biochem_fan
      License: dual license of MIT or LGPL3

  Contributors:
    Robert Hanson for parseXYZ, deferred instantiation

  This program uses
      Three.js
         https://github.com/mrdoob/three.js
         Copyright (c) 2010-2012 three.js Authors. All rights reserved.
      jQuery
         http://jquery.org/
         Copyright (c) 2011 John Resig
 */



(function (window, undefined) {
    'use strict';
    var hasCanvas, hasWebgl, TV3, TF3, TCo, isNotSolvent;

    hasCanvas = (function () {
        var canvas  = document.createElement('canvas');
        return !!(canvas.getContext && canvas.getContext('2d'));
    }());

    hasWebgl = (function () {
        if (!hasCanvas) { return false; }
        var canvas  = document.createElement('canvas');
        try {
            return !!(window.WebGLRenderingContext && (canvas.getContext('webgl') || canvas.getContext('experimental-webgl')));
        } catch (e) {
            return false;
        }
    }());

    // Workaround for Intel GMA series (gl_FrontFacing causes compilation error)
    THREE.ShaderLib.lambert.fragmentShader = THREE.ShaderLib.lambert.fragmentShader.replace("gl_FrontFacing", "true");
    THREE.ShaderLib.lambert.vertexShader = THREE.ShaderLib.lambert.vertexShader.replace(/\}$/, "#ifdef DOUBLE_SIDED\n if (transformedNormal.z < 0.0) vLightFront = vLightBack;\n #endif\n }");

    TV3 = THREE.Vector3;
    TF3 = THREE.Face3;
    TCo = THREE.Color;

    THREE.Geometry.prototype.colorAll = function (color) {
        this.faces.forEach(function (face) { face.color = color; });
    };

    THREE.Matrix4.prototype.isIdentity = function() {
        var i, j;
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) { 
                if (this.elements[i * 4 + j] != (i === j) ? 1 : 0) { return false; }
            }
        }
        return true;
    };



    function GLmol(queryselector, suppressAutoload, fallback) {
        if (!queryselector) {
            return false;
        }
        this.queryselector = '#' + queryselector;
        this.aaScale = 1; // or 2

        this.container = $(this.queryselector);
        this.WIDTH = this.container.width() * this.aaScale;
        this.HEIGHT = this.container.height() * this.aaScale;
        this.ASPECT = this.WIDTH / this.HEIGHT;

        if (hasWebgl) {
            this.renderer = new THREE.WebGLRenderer({antialias: true});
        } else if (hasCanvas && fallback) {
            this.renderer = new THREE.CanvasRenderer();
        } else {
            throw new Error("no suitable renderer");
        }

        this.renderer.sortObjects = false; // hopefully improve performance
        // 'antialias: true' now works in Firefox too!
        // setting this.aaScale = 2 will enable antialias in older Firefox but GPU load increases.
        this.renderer.domElement.style.width = "100%";
        this.renderer.domElement.style.height = "100%";
        this.container.append(this.renderer.domElement);
        this.renderer.setSize(this.WIDTH, this.HEIGHT);

        this.camera = new THREE.PerspectiveCamera(20, this.ASPECT, 1, 800); // will be updated anyway
        this.camera.position = new TV3(0, 0, this.CAMERA_Z);
        this.camera.lookAt(new TV3(0, 0, 0));
        this.perspectiveCamera = this.camera;
        this.orthoscopicCamera = new THREE.OrthographicCamera();
        this.orthoscopicCamera.position.z = this.CAMERA_Z;
        this.orthoscopicCamera.lookAt(new TV3(0, 0, 0));

        this.attachResizeEvent();
        this.onWindowResize();

        // UI variables
        this.cq = new THREE.Quaternion(1, 0, 0, 0);
        this.dq = new THREE.Quaternion(1, 0, 0, 0);
        this.isDragging = false;
        this.mouseStartX = 0;
        this.mouseStartY = 0;
        this.currentModelPos = 0;
        this.cz = 0;

        this.enableMouse(); // TODO: think about using default THREE.js controls

        if (!suppressAutoload) {
            this.loadMolecule();
        }

    };

    GLmol.prototype.Nucleotides = ['  G', '  A', '  T', '  C', '  U', ' DG', ' DA', ' DT', ' DC', ' DU'];
    GLmol.prototype.ElementColors = {
        "H": 0xCCCCCC,
        "C": 0xAAAAAA,
        "O": 0xCC0000,
        "N": 0x0000CC,
        "S": 0xCCCC00,
        "P": 0x6622CC,
        "F": 0x00CC00,
        "Cl": 0x00CC00,
        "Br": 0x882200,
        "I": 0x6600AA,
        "Fe": 0xCC6600,
        "Ca": 0x8888AA
    };
    // Reference: A. Bondi, J. Phys. Chem., 1964, 68, 441.
    GLmol.prototype.vdwRadii = {
        "H": 1.2,
        "He": 1.4,
        "Li": 1.82,
        "Be": 1.53,
        "B": 1.92,
        "C": 1.7,
        "N": 1.55,
        "O": 1.52,
        "F": 1.47,
        "Ne": 1.54,
        "Na": 2.27,
        "Mg": 1.73,
        "Al": 1.84,
        "Si": 2.10,
        "P": 1.80,
        "S": 1.80,
        "Cl": 1.75,
        "Ar": 1.88,
        "K": 2.75,
        "Ca": 2.31,
        "Sc": 2.11,
        // Ti
        // V
        // Cr
        // Mn
        // Fe
        // Co
        "Ni": 1.63,
        "Cu": 1.4,
        "Zn": 1.39,
        "Ga": 1.87,
        "Ge": 2.11,
        "As": 1.85,
        "Se": 1.90,
        "Br": 1.85,
        "Kr": 2.02,
        "Rb": 3.03,
        "Sr": 2.49,
        // Y,
        // Zr
        // Nb
        // Mo
        // Tc
        // Ru
        // Rh
        "Pd": 1.63,
        "Ag": 1.72,
        "Cd": 1.58,
        "In": 1.93,
        "Sn": 2.17,
        "Sb": 2.06,
        "Te": 2.06,
        "I": 1.98,
        "Xe": 2.16,
        "Cs": 3.43,
        "Ba": 2.68,
        // La
        // Ce
        // Pr
        // Nd
        // Pm
        // Sm
        // Eu
        // Gd
        // Tb
        // Dy
        // Ho
        // Er
        // Tm
        // Yb
        // Lu
        // Hf
        // Ta
        // W
        // Re
        // Os
        // Ir
        "Pt": 1.75,
        "Au": 1.66,
        "Hg": 1.55,
        "Tl": 1.96,
        "Pb": 2.02,
        "Bi": 2.07,
        "Po": 1.97,
        "At": 2.02,
        "Rn": 2.20,
        "Fr": 3.48,
        "Ra": 2.83,
        //Ac
        //Th
        //Pa
        "U": 1.86,
        //Np
        //Pu
        //Am
        //Cm
        //Bk
        //Cf
        //Es
        //Fm
        //Md
        //No
        //Lr
        //Rf
        //Db
        //Sg
        //Bh
        //Hs
        //Mt
        //Ds
        //Rg
        //Cn
        //Uut
        //Fl
        //Uup
        //Lv
        //Uus
        //Uuo
    };

    GLmol.prototype.NEAR = 1;
    GLmol.prototype.FAR = 800;
    GLmol.prototype.CAMERA_Z = -150;
    GLmol.prototype.aaScale = 1;
    GLmol.prototype.scene = null;
    GLmol.prototype.rotationGroup = null; // which contains modelGroup
    GLmol.prototype.labelGroup = null; // which contains modelGroup
    GLmol.prototype.modelGroup = null;
    GLmol.prototype.bgColor = 0x000000;
    GLmol.prototype.fov = 20;
    GLmol.prototype.fogStart = 0.4;
    GLmol.prototype.slabNear = -50; // relative to the center of rotationGroup
    GLmol.prototype.slabFar = +50;
    // Default values
    GLmol.prototype.sphereRadius = 1.5;
    GLmol.prototype.cylinderRadius = 0.4;
    GLmol.prototype.lineWidth = 1.5 * GLmol.prototype.aaScale;
    GLmol.prototype.curveWidth = 3 * GLmol.prototype.aaScale;
    GLmol.prototype.defaultColor = 0xCCCCCC;
    GLmol.prototype.sphereQuality = 16; //16;
    GLmol.prototype.cylinderQuality = 16; //8;
    GLmol.prototype.axisDIV = 5; // 3 still gives acceptable quality
    GLmol.prototype.strandDIV = 6;
    GLmol.prototype.nucleicAcidStrandDIV = 4;
    GLmol.prototype.tubeDIV = 8;
    GLmol.prototype.coilWidth = 0.3;
    GLmol.prototype.helixSheetWidth = 1.3;
    GLmol.prototype.nucleicAcidWidth = 0.8;
    GLmol.prototype.thickness = 0.4;
    GLmol.prototype.protein = {sheet: [], helix: [], biomtChains: '', biomtMatrices: [], symMat: [], pdbID: '', title: ''};
    GLmol.prototype.atoms = [];

    GLmol.prototype.onWindowResize = function (e) {
        this.WIDTH = this.container.width() * this.aaScale;
        this.HEIGHT = this.container.height() * this.aaScale;
        this.ASPECT = this.WIDTH / this.HEIGHT;
        this.renderer.setSize(this.WIDTH, this.HEIGHT);
        this.camera.aspect = this.ASPECT;
        this.camera.updateProjectionMatrix();
        this.show();
    }
    GLmol.prototype.attachResizeEvent = function () {
        $(window).resize(function (e) { // only window can capture resize event
            this.onWindowResize(e); // needed in fullscreen mode too
        }.bind(this));
    };

    GLmol.prototype.setupLights = function (scene) {
        var directionalLight, ambientLight;

        directionalLight =  new THREE.DirectionalLight(0xFFFFFF); //MAGIC
        directionalLight.position = new TV3(0.2, 0.2, -1).normalize();
        directionalLight.intensity = 1.2;
        scene.add(directionalLight);

        ambientLight = new THREE.AmbientLight(0x202020); //MAGIC
        scene.add(ambientLight);
    };

    GLmol.prototype.parseSDF = function (str) {
        var atoms = this.atoms,
            protein = this.protein,
            lines = str.split("\n"),
            offset = 4,
            atomCount,
            bondCount,
            i,
            line,
            atom,
            from,
            to,
            order;

        if (lines.length < 4) { return; }
        atomCount = parseInt(lines[3].substr(0, 3), 10);
        if (isNaN(atomCount) || atomCount <= 0) { return; }
        bondCount = parseInt(lines[3].substr(3, 3), 10);
        if (lines.length < 4 + atomCount + bondCount) { return; }

        for (i = 1; i <= atomCount; i++) {
            line = lines[offset];
            atom = {};
            offset++;
            atom.serial = i;
            atom.x = parseFloat(line.substr(0, 10));
            atom.y = parseFloat(line.substr(10, 10));
            atom.z = parseFloat(line.substr(20, 10));
            atom.hetflag = true;
            atom.atom = atom.elem = line.substr(31, 3).replace(/ /g, "");
            atom.bonds = [];
            atom.bondOrder = [];
            atoms[i] = atom;
        }

        for (i = 1; i <= bondCount; i++) {
            line = lines[offset];
            from = parseInt(line.substr(0, 3), 10);
            to = parseInt(line.substr(3, 3), 10);
            order = parseInt(line.substr(6, 3), 10);
            offset++;
            atoms[from].bonds.push(to);
            atoms[from].bondOrder.push(order);
            atoms[to].bonds.push(from);
            atoms[to].bondOrder.push(order);
        }

        protein.smallMolecule = true;
        return true;
    };

    GLmol.prototype.parseXYZ = function (str) {
        var atoms = this.atoms,
            protein = this.protein,
            lines = str.split("\n"),
            atomCount,
            offset,
            i,
            j,
            line,
            tokens,
            atom;

        if (lines.length < 3) { return; }
        atomCount = parseInt(lines[0].substr(0, 3), 10);
        if (isNaN(atomCount) || atomCount <= 0) { return; }
        if (lines.length < atomCount + 2) { return; }
        offset = 2;
        for (i = 1; i <= atomCount; i++) {
            line = lines[offset++];
            tokens = line.replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
            console.log(tokens);
            atom = {};
            atom.serial = i;
            atom.atom = atom.elem = tokens[0];
            atom.x = parseFloat(tokens[1]);
            atom.y = parseFloat(tokens[2]);
            atom.z = parseFloat(tokens[3]);
            atom.hetflag = true;
            atom.bonds = [];
            atom.bondOrder = [];
            atoms[i] = atom;
        }
        for (i = 1; i < atomCount; i++) { // hopefully XYZ is small enough
            for (j = i + 1; j <= atomCount; j++) {
                if (this.isConnected(atoms[i], atoms[j])) {
                    atoms[i].bonds.push(j);
                    atoms[i].bondOrder.push(1);
                    atoms[j].bonds.push(i);
                    atoms[j].bondOrder.push(1);
                }
            }
        }
        protein.smallMolecule = true;
        return true;
    };

    GLmol.prototype.parsePDB2 = function (str) {
        var atoms = this.atoms,
            protein = this.protein,
            molID,
            atoms_cnt = 0,
            lines = str.split("\n"),
            i,
            j,
            line,
            recordName,
            atom,
            resn,
            chain,
            resi,
            x,
            y,
            z,
            hetflag,
            elem,
            serial,
            altLoc,
            b,
            startChain,
            startResi,
            endChain,
            endResi,
            from,
            to,
            type,
            n,
            m,
            found;


        for (i = 0; i < lines.length; i++) {
            line = lines[i].replace(/^\s*/, ''); // remove indent
            recordName = line.substr(0, 6);
            if (recordName === 'ATOM  ' || recordName === 'HETATM') {
                altLoc = line.substr(16, 1);
                if (altLoc !== ' ' && altLoc !== 'A') { continue; } // FIXME: ad hoc
                serial = parseInt(line.substr(6, 5), 10);
                atom = line.substr(12, 4).replace(/ /g, "");
                resn = line.substr(17, 3);
                chain = line.substr(21, 1);
                resi = parseInt(line.substr(22, 5), 10);
                x = parseFloat(line.substr(30, 8));
                y = parseFloat(line.substr(38, 8));
                z = parseFloat(line.substr(46, 8));
                b = parseFloat(line.substr(60, 8));
                elem = line.substr(76, 2).replace(/ /g, "");
                if (elem === '') { // for some incorrect PDB files
                    elem = line.substr(12, 4).replace(/ /g, "");
                }

                hetflag = line[0] === 'H';

                atoms[serial] = {'resn': resn, 
                                 'x': x,
                                 'y': y,
                                 'z': z,
                                 'elem': elem,
                                 'hetflag': hetflag, 
                                 'chain': chain, 
                                 'resi': resi, 
                                 'serial': serial, 
                                 'atom': atom,
                                 'bonds': [], 
                                 'ss': 'c',
                                 'color': 0xFFFFFF,
                                 'bondOrder': [],
                                 'b': b
                               /*', altLoc': altLoc*/};
            } else if (recordName === 'SHEET ') {
                startChain = line.substr(21, 1);
                startResi = parseInt(line.substr(22, 4), 10);
                endChain = line.substr(32, 1);
                endResi = parseInt(line.substr(33, 4), 10);
                protein.sheet.push([startChain, startResi, endChain, endResi]);
            } else if (recordName === 'CONECT') {
    // MEMO: We don't have to parse SSBOND, LINK because both are also
    // described in CONECT. But what about 2JYT???
                from = parseInt(line.substr(6, 5), 10);
                for (j = 0; j < 4; j++) {
                    to = parseInt(line.substr([11, 16, 21, 26][j], 5), 10);
                    if (isNaN(to)) { continue; }
                    if (atoms[from] !== undefined) {
                        atoms[from].bonds.push(to);
                        atoms[from].bondOrder.push(1);
                    }
                }
            } else if (recordName === 'HELIX ') {
                startChain = line.substr(19, 1);
                startResi = parseInt(line.substr(21, 4), 10);
                endChain = line.substr(31, 1);
                endResi = parseInt(line.substr(33, 4), 10);
                protein.helix.push([startChain, startResi, endChain, endResi]);
            } else if (recordName === 'CRYST1') {
                protein.a = parseFloat(line.substr(6, 9));
                protein.b = parseFloat(line.substr(15, 9));
                protein.c = parseFloat(line.substr(24, 9));
                protein.alpha = parseFloat(line.substr(33, 7));
                protein.beta = parseFloat(line.substr(40, 7));
                protein.gamma = parseFloat(line.substr(47, 7));
                protein.spacegroup = line.substr(55, 11);
                this.defineCell();
            } else if (recordName === 'REMARK') {
                type = parseInt(line.substr(7, 3), 10);
                if (type === 290 && line.substr(13, 5) === 'SMTRY') {
                    n = parseInt(line[18], 10) - 1;
                    m = parseInt(line.substr(21, 2), 10);
                    if (!protein.symMat[m]) {
                        protein.symMat[m] = new THREE.Matrix4().identity();
                    }
                    protein.symMat[m].elements[n] = parseFloat(line.substr(24, 9));
                    protein.symMat[m].elements[n + 4] = parseFloat(line.substr(34, 9));
                    protein.symMat[m].elements[n + 8] = parseFloat(line.substr(44, 9));
                    protein.symMat[m].elements[n + 12] = parseFloat(line.substr(54, 10));
                } else if (type === 350 && line.substr(13, 5) === 'BIOMT') {
                    n = parseInt(line[18], 10) - 1;
                    m = parseInt(line.substr(21, 2), 10);
                    if (!protein.biomtMatrices[m]) {
                        protein.biomtMatrices[m] = new THREE.Matrix4().identity();
                    }
                    protein.biomtMatrices[m].elements[n] = parseFloat(line.substr(24, 9));
                    protein.biomtMatrices[m].elements[n + 4] = parseFloat(line.substr(34, 9));
                    protein.biomtMatrices[m].elements[n + 8] = parseFloat(line.substr(44, 9));
                    protein.biomtMatrices[m].elements[n + 12] = parseFloat(line.substr(54, 10));
                } else if (type === 350 && line.substr(11, 11) === 'BIOMOLECULE') {
                    protein.biomtMatrices = [];
                    protein.biomtChains = '';
                } else if (type === 350 && line.substr(34, 6) === 'CHAINS') {
                    protein.biomtChains += line.substr(41, 40);
                }
            } else if (recordName === 'HEADER') {
                protein.pdbID = line.substr(62, 4);
            } else if (recordName === 'TITLE ') {
                if (!protein.title) {
                    protein.title = "";
                } else {
                    protein.title += line.substr(10, 70) + "\n"; // CHECK: why 60 is not enough???
                }
            } else if (recordName === 'COMPND') {
                console.warn("CMPND record unimplemented"); // TODO: Implement me!
            }
        }

       // Assign secondary structures
        for (i = 0; i < atoms.length; i++) {
            atom = atoms[i];
            if (!atom) { continue; }

            found = false;
            // MEMO: Can start chain and end chain differ?
            for (j = 0; j < protein.sheet.length; j++) {
                if (atom.chain !== protein.sheet[j][0]) { continue; }
                if (atom.resi < protein.sheet[j][1]) { continue; }
                if (atom.resi > protein.sheet[j][3]) { continue; }
                atom.ss = 's';
                if (atom.resi === protein.sheet[j][1]) { atom.ssbegin = true; }
                if (atom.resi === protein.sheet[j][3]) { atom.ssend = true; }
            }
            for (j = 0; j < protein.helix.length; j++) {
                if (atom.chain !== protein.helix[j][0]) { continue; }
                if (atom.resi < protein.helix[j][1]) { continue; }
                if (atom.resi > protein.helix[j][3]) { continue; }
                atom.ss = 'h';
                if (atom.resi === protein.helix[j][1]) {
                    atom.ssbegin = true;
                } else if (atom.resi === protein.helix[j][3]) {
                    atom.ssend = true;
                }
            }
        }
        protein.smallMolecule = false;
        return true;
    };

    // Catmull-Rom subdivision
    GLmol.prototype.subdivide = function (_points, DIV) { // points as Vector3
        var ret = [],
            points = _points,
            i,
            lim,
            p0,
            p1,
            p2,
            p3,
            size,
            v0,
            v1,
            t,
            x,
            y,
            z,
            j;

        points = []; // Smoothing test
        points.push(_points[0]);
        for (i = 1, lim = _points.length - 1; i < lim; i++) {
            p1 = _points[i];
            p2 = _points[i + 1];
            if (p1.smoothen) {
                points.push(new TV3((p1.x + p2.x) / 2, (p1.y + p2.y) / 2, (p1.z + p2.z) / 2));
            } else {
                points.push(p1);
            }
        }
        points.push(_points[_points.length - 1]);

        for (i = -1, size = points.length; i <= size - 3; i++) {
            p0 = points[(i === -1) ? 0 : i];
            p1 = points[i + 1];
            p2 = points[i + 2];
            p3 = points[(i === size - 3) ? size - 1 : i + 3];
            v0 = new TV3().sub(p2, p0).multiplyScalar(0.5);
            v1 = new TV3().sub(p3, p1).multiplyScalar(0.5);
            for (j = 0; j < DIV; j++) {
                t = 1.0 / DIV * j;
                x = p1.x + t * v0.x
                     + t * t * (-3 * p1.x + 3 * p2.x - 2 * v0.x - v1.x)
                     + t * t * t * (2 * p1.x - 2 * p2.x + v0.x + v1.x);
                y = p1.y + t * v0.y
                     + t * t * (-3 * p1.y + 3 * p2.y - 2 * v0.y - v1.y)
                     + t * t * t * (2 * p1.y - 2 * p2.y + v0.y + v1.y);
                z = p1.z + t * v0.z
                     + t * t * (-3 * p1.z + 3 * p2.z - 2 * v0.z - v1.z)
                     + t * t * t * (2 * p1.z - 2 * p2.z + v0.z + v1.z);
                ret.push(new TV3(x, y, z));
            }
        }
        ret.push(points[points.length - 1]);
        return ret;
    };

    GLmol.prototype.getLambertMesh = function (geometry, options) {
        var material = new THREE.MeshLambertMaterial({
            vertexColors: THREE.FaceColors,
            side: THREE.DoubleSide
        });
        return this.getMesh(material, geometry, options);
    };


    GLmol.prototype.getMesh = function (material, geometry, options) {
        material.setValues(options);
        return new THREE.Mesh(geometry, material);
    };

    GLmol.prototype.drawAtoms = function (geometry, group, atomlist, defaultRadius, forceDefault, scale) {
        atomlist.filter(isNotUndefined).forEach(function (atom) {
            var sphere,
                r = defaultRadius;

            if (!forceDefault) {
                scale = scale || 1;
                r = (this.vdwRadii[atom.elem] || defaultRadius) * scale;
            }

            sphere = this.getLambertMesh(geometry, {color: atom.color});
            sphere.scale.x = sphere.scale.y = sphere.scale.z = r;
            sphere.position.x = atom.x;
            sphere.position.y = atom.y;
            sphere.position.z = atom.z;

            group.add(sphere);
        }.bind(this));
    };
    GLmol.prototype.drawAtomsAsSphere = function (group, atomlist, defaultRadius, forceDefault, scale) {
        var geometry = new THREE.SphereGeometry(1, this.sphereQuality, this.sphereQuality); // r, seg, ring
        this.drawAtoms(geometry, group, atomlist, defaultRadius, forceDefault, scale);
    };

    // about two times faster than sphere when div = 2
    GLmol.prototype.drawAtomsAsIcosahedron = function (group, atomlist, defaultRadius, forceDefault, scale) {
        var geometry = this.IcosahedronGeometry();
        this.drawAtoms(geometry, group, atomlist, defaultRadius, forceDefault, scale);
    };

    GLmol.prototype.isConnected = function (atom1, atom2) {
        var s = atom1.bonds.indexOf(atom2.serial),
            distSquared;
        if (s !== -1) { return atom1.bondOrder[s]; }

        if (this.protein.smallMolecule && (atom1.hetflag || atom2.hetflag)) { return 0; } // CHECK: or should I ?

        distSquared = (atom1.x - atom2.x) * (atom1.x - atom2.x) +
                      (atom1.y - atom2.y) * (atom1.y - atom2.y) +
                      (atom1.z - atom2.z) * (atom1.z - atom2.z);

    //   if (atom1.altLoc != atom2.altLoc) return 0;
        if (isNaN(distSquared)) { return 0; }
        if (distSquared < 0.5) { return 0; } // maybe duplicate position.

        if (distSquared > 1.3 && (atom1.elem === 'H' || atom2.elem === 'H' || atom1.elem === 'D' || atom2.elem === 'D')) { return 0; }
        if (distSquared < 3.42 && (atom1.elem === 'S' || atom2.elem === 'S')) { return 1; }
        if (distSquared > 2.78) { return 0; }
        return 1;
    };

    GLmol.prototype.drawBondAsStickSub = function (group, atom1, atom2, bondR, order) {
        var delta, tmp, p1, p2, mp, mp3, c1, c2;

        p1 = new TV3(atom1.x, atom1.y, atom1.z);
        p2 = new TV3(atom2.x, atom2.y, atom2.z);
        mp = p1.clone().addSelf(p2).multiplyScalar(0.5);

        c1 = new TCo(atom1.color);
        c2 = new TCo(atom2.color);
        if (order === 1 || order === 3) {
            this.drawCylinder(group, p1, mp, bondR, atom1.color);
            this.drawCylinder(group, p2, mp, bondR, atom2.color);
        }
        if (order > 1) {
            delta = this.calcBondDelta(atom1, atom2, bondR * 2.3);
            tmp = mp.clone().addSelf(delta);
            this.drawCylinder(group, p1.clone().addSelf(delta), tmp, bondR, atom1.color);
            this.drawCylinder(group, p2.clone().addSelf(delta), tmp, bondR, atom2.color);
            tmp = mp.clone().subSelf(delta);
            this.drawCylinder(group, p1.clone().subSelf(delta), tmp, bondR, atom1.color);
            this.drawCylinder(group, p2.clone().subSelf(delta), tmp, bondR, atom2.color);
        }
    };

    GLmol.prototype.drawBondsAsStick = function (group, atomlist, bondR, atomR, ignoreNonbonded, multipleBonds, scale) {
        var nAtoms = atomlist.length,
            mp,
            forSpheres = [],
            _i,
            i,
            _j,
            j,
            atom1,
            atom2,
            order;

        if (multipleBonds) {
            bondR /= 2.5;
        }

        for (_i = 0; _i < nAtoms; _i++) {
            atom1 = atomlist[_i];
            if (!atom1) { continue; }
            for (_j = _i + 1; _j < _i + 30 && _j < nAtoms; _j++) {
                atom2 = atomlist[_j];
                if (!atom2) {
                    continue;
                }
                order = this.isConnected(atom1, atom2);
                if (order === 0) { continue; }
                atom1.connected = atom2.connected = true;
                this.drawBondAsStickSub(group, atom1, atom2, bondR, (!!multipleBonds) ? order : 1);
            }
            for (_j = 0; _j < atom1.bonds.length; _j++) {
                j = atom1.bonds[_j];
                if (j < i + 30) { continue; } // be conservative!
                if (atomlist.indexOf(j) === -1) { continue; }
                atom2 = this.atoms[j];
                if (!atom2) { continue; }
                atom1.connected = atom2.connected = true;
                this.drawBondAsStickSub(group, atom1, atom2, bondR, (!!multipleBonds) ? atom1.bondOrder[_j] : 1);
            }

            if (atom1.connected) {
                forSpheres.push(atom1);
            }
        }
        this.drawAtomsAsSphere(group, forSpheres, atomR, !scale, scale);
    };

    GLmol.prototype.defineCell = function () {
        var p = this.protein;
        if (!p.a) { return; }

        p.ax = p.a;
        p.ay = 0;
        p.az = 0;
        p.bx = p.b * Math.cos(Math.PI / 180.0 * p.gamma);
        p.by = p.b * Math.sin(Math.PI / 180.0 * p.gamma);
        p.bz = 0;
        p.cx = p.c * Math.cos(Math.PI / 180.0 * p.beta);
        p.cy = p.c * (Math.cos(Math.PI / 180.0 * p.alpha) -
                   Math.cos(Math.PI / 180.0 * p.gamma)
                 * Math.cos(Math.PI / 180.0 * p.beta)
                 / Math.sin(Math.PI / 180.0 * p.gamma));
        p.cz = Math.sqrt(p.c * p.c * Math.sin(Math.PI / 180.0 * p.beta)
                   * Math.sin(Math.PI / 180.0 * p.beta) - p.cy * p.cy);
    };

    GLmol.prototype.drawUnitcell = function (group) {
        var p = this.protein,
            vertices,
            edges,
            geo,
            lineMaterial,
            line,
            i;

        if (!p.a) { return; }

        vertices = [[0, 0, 0],
                    [p.ax, p.ay, p.az],
                    [p.bx, p.by, p.bz],
                    [p.ax + p.bx, p.ay + p.by, p.az + p.bz],
                    [p.cx, p.cy, p.cz],
                    [p.cx + p.ax, p.cy + p.ay,  p.cz + p.az],
                    [p.cx + p.bx, p.cy + p.by, p.cz + p.bz],
                    [p.cx + p.ax + p.bx, p.cy + p.ay + p.by, p.cz + p.az + p.bz]];

        edges = [0, 1, 0, 2, 1, 3, 2, 3, 4, 5, 4, 6, 5, 7, 6, 7, 0, 4, 1, 5, 2, 6, 3, 7];

        geo = new THREE.Geometry();
        geo.vertices = edges.map(function (edge) {
            return (new TV3(vertices[edge][0], vertices[edge][1], vertices[edge][2]));
        });

        lineMaterial = new THREE.LineBasicMaterial({linewidth: 1, color: 0xcccccc}); // MAGIC
        line = new THREE.Line(geo, lineMaterial);
        line.type = THREE.LinePieces;
        group.add(line);
    };

    // TODO: Find inner side of a ring
    GLmol.prototype.calcBondDelta = function (atom1, atom2, sep) {
        var dot,
            axis = new TV3(atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z).normalize(),
            found = null,
            i,
            tmp,
            atom,
            delta;

        for (i = 0; i < atom1.bonds.length && !found; i++) {
            atom = this.atoms[atom1.bonds[i]];
            if (!atom) {
                continue;
            }
            if (atom.serial !== atom2.serial && atom.elem !== 'H') {
                found = atom;
            }
        }
        for (i = 0; i < atom2.bonds.length && !found; i++) {
            atom = this.atoms[atom2.bonds[i]];
            if (!atom) {
                continue;
            }
            if (atom.serial !== atom1.serial && atom.elem !== 'H') {
                found = atom;
            }
        }
        if (found) {
            tmp = new TV3(atom1.x - found.x, atom1.y - found.y, atom1.z - found.z).normalize();
            dot = tmp.dot(axis);
            delta = new TV3(tmp.x - axis.x * dot, tmp.y - axis.y * dot, tmp.z - axis.z * dot);
        }
        if (!found || Math.abs(dot - 1) < 0.001 || Math.abs(dot + 1) < 0.001) {
            if (axis.x < 0.01 && axis.y < 0.01) {
                delta = new TV3(0, -axis.z, axis.y);
            } else {
                delta = new TV3(-axis.y, axis.x, 0);
            }
        }
        delta.normalize().multiplyScalar(sep);
        return delta;
    };

    GLmol.prototype.drawBondsAsLineSub = function (geo, atom1, atom2, order) {
        var delta,
            tmp,
            vs = geo.vertices,
            cs = geo.colors,
            p1,
            p2,
            mp,
            c1,
            c2;

        if (order > 1) { delta = this.calcBondDelta(atom1, atom2, 0.15); }
        p1 = new TV3(atom1.x, atom1.y, atom1.z);
        p2 = new TV3(atom2.x, atom2.y, atom2.z);
        mp = p1.clone().addSelf(p2).multiplyScalar(0.5);

        c1 = new TCo(atom1.color);
        c2 = new TCo(atom2.color);
        if (order === 1 || order === 3) {
            geo.vertices.push(p1);
            geo.colors.push(c1);
            geo.vertices.push(mp);
            geo.colors.push(c1);
            geo.vertices.push(p2);
            geo.colors.push(c2);
            geo.vertices.push(mp);
            geo.colors.push(c2);
        }
        if (order > 1) {
            geo.vertices.push(p1.clone().addSelf(delta));
            geo.colors.push(c1);
            geo.vertices.push(tmp = mp.clone().addSelf(delta));
            geo.colors.push(c1);
            geo.vertices.push(p2.clone().addSelf(delta));
            geo.colors.push(c2);
            geo.vertices.push(tmp);
            geo.colors.push(c2);
            geo.vertices.push(p1.clone().subSelf(delta));
            geo.colors.push(c1);
            geo.vertices.push(tmp = mp.clone().subSelf(delta));
            geo.colors.push(c1);
            geo.vertices.push(p2.clone().subSelf(delta));
            geo.colors.push(c2);
            geo.vertices.push(tmp);
            geo.colors.push(c2);
        }
    };

    GLmol.prototype.drawBondsAsLine = function (group, atomlist, lineWidth) {
        var geo = new THREE.Geometry(),
            nAtoms = atomlist.length,
            _i,
            i,
            _j,
            j,
            atom1,
            atom2,
            order,
            lineMaterial,
            line;


        for (_i = 0; _i < nAtoms; _i++) {
            atom1 = atomlist[_i];
            if (!atom1) { continue; }
            for (_j = _i + 1; _j < _i + 30 && _j < nAtoms; _j++) {
                atom2 = atomlist[_j];
                if (!atom2) { continue; }
                order = this.isConnected(atom1, atom2);
                if (order === 0) { continue; }

                this.drawBondsAsLineSub(geo, atom1, atom2, order);
            }
            for (_j = 0; _j < atom1.bonds.length; _j++) {
                j = atom1.bonds[_j];
                if (j < i + 30) { continue; } // be conservative!
                if (atomlist.indexOf(j) === -1) { continue; }
                atom2 = this.atoms[j];
                if (!atom2) { continue; }
                this.drawBondsAsLineSub(geo, atom1, atom2, atom1.bondOrder[_j]);
            }
        }

        lineMaterial = new THREE.LineBasicMaterial({linewidth: lineWidth});
        lineMaterial.vertexColors = true;

        line = new THREE.Line(geo, lineMaterial);
        line.type = THREE.LinePieces;
        group.add(line);
    };

    GLmol.prototype.drawSmoothCurve = function (group, _points, width, colors, div) {
        if (!_points || _points.length === 0) { return; }

        div = (div) ? 5 : div;

        var geo = new THREE.Geometry(),
            points = this.subdivide(_points, div),
            i,
            lineMaterial,
            line;

        for (i = 0; i < points.length; i++) {
            geo.vertices.push(points[i]);
            geo.colors.push(new TCo(colors[(i === 0) ? 0 : Math.round((i - 1) / div)]));
        }

        //line = this.getLineMesh({linewidth: width, type: THREE.LineStrip});
        lineMaterial = new THREE.LineBasicMaterial({linewidth: width});
        lineMaterial.vertexColors = true;
        line = new THREE.Line(geo, lineMaterial);
        line.type = THREE.LineStrip;
        group.add(line);
    };

    GLmol.prototype.drawAsCross = function (group, atomlist, delta) {
        var geo = new THREE.Geometry(),
            points = [[delta, 0, 0], [-delta, 0, 0], [0, delta, 0], [0, -delta, 0], [0, 0, delta], [0, 0, -delta]],
            c,
            lineMaterial,
            line,
            i,
            lim,
            atom,
            j;

        atomlist = atomlist.filter(isNotUndefined);

        for (i = 0, lim = atomlist.length; i < lim; i++) {
            atom = atomlist[i];

            c = new TCo(atom.color);
            for (j = 0; j < 6; j++) {
                geo.vertices.push(new TV3(atom.x + points[j][0], atom.y + points[j][1], atom.z + points[j][2]));
                geo.colors.push(c);
            }
        }

        lineMaterial = new THREE.LineBasicMaterial({linewidth: this.lineWidth});
        lineMaterial.vertexColors = true;
        line = new THREE.Line(geo, lineMaterial, THREE.LinePieces);
        group.add(line);
    };

    // FIXME: Winkled...
    GLmol.prototype.drawSmoothTube = function (group, _points, colors, radii) {
        if (_points.length < 2) { return; }

        var circleDiv = this.tubeDIV,
            axisDiv = this.axisDIV,
            geo = new THREE.Geometry(),
            points = this.subdivide(_points, axisDiv),
            prevAxis1 = new TV3(),
            prevAxis2,
            r,
            idx,
            floored,
            tmp,
            delta,
            axis1,
            axis2,
            i,
            lim,
            j,
            angle,
            c,
            s,
            offset,
            reg,
            r1,
            r2,
            mat,
            mesh;

        for (i = 0, lim = points.length; i < lim; i++) {
            idx = (i - 1) / axisDiv;
            if (i === 0) {
                r = radii[0];
            } else {
                if (idx % 1 === 0) {
                    r = radii[idx];
                } else {
                    floored = Math.floor(idx);
                    tmp = idx - floored;
                    r = radii[floored] * tmp + radii[floored + 1] * (1 - tmp);
                }
            }

            if (i < lim - 1) {
                delta = new TV3().sub(points[i], points[i + 1]);
                axis1 = new TV3(0, -delta.z, delta.y).normalize().multiplyScalar(r);
                axis2 = new TV3().cross(delta, axis1).normalize().multiplyScalar(r);
    //          var dir = 1, offset = 0;
                if (prevAxis1.dot(axis1) < 0) {
                    axis1.negate();
                    axis2.negate();  //dir = -1;//offset = 2 * Math.PI / axisDiv;
                }
                prevAxis1 = axis1;
                prevAxis2 = axis2;
            } else {
                axis1 = prevAxis1;
                axis2 = prevAxis2;
            }

            for (j = 0; j < circleDiv; j++) {
                angle = 2 * Math.PI / circleDiv * j; //* dir  + offset;
                c = Math.cos(angle);
                s = Math.sin(angle);
                geo.vertices.push(new TV3(
                    points[i].x + c * axis1.x + s * axis2.x,
                    points[i].y + c * axis1.y + s * axis2.y,
                    points[i].z + c * axis1.z + s * axis2.z
                ));
            }
        }

        offset = 0;
        for (i = 0, lim = points.length - 1; i < lim; i++) {
            c =  new TCo(colors[Math.round((i - 1) / axisDiv)]);

            reg = 0;
            r1 = new TV3().sub(geo.vertices[offset], geo.vertices[offset + circleDiv]).lengthSq();
            r2 = new TV3().sub(geo.vertices[offset], geo.vertices[offset + circleDiv + 1]).lengthSq();
            if (r1 > r2) {r1 = r2; reg = 1; }
            for (j = 0; j < circleDiv; j++) {
                geo.faces.push(new TF3(offset + j, offset + (j + reg) % circleDiv + circleDiv, offset + (j + 1) % circleDiv));
                geo.faces.push(new TF3(offset + (j + 1) % circleDiv, offset + (j + reg) % circleDiv + circleDiv, offset + (j + reg + 1) % circleDiv + circleDiv));
                geo.faces[geo.faces.length - 2].color = c;
                geo.faces[geo.faces.length - 1].color = c;
            }
            offset += circleDiv;
        }
        geo.computeFaceNormals();
        geo.computeVertexNormals(false);

        mesh = this.getLambertMesh(geo);
        group.add(mesh);
    };


    GLmol.prototype.drawMainchainCurve = function (group, atomlist, curveWidth, atomName, div) {
        var points = [],
            colors = [],
            currentChain,
            currentResi,
            i,
            atom;

        div = div || 5;

        atomlist = atomlist.filter(isNotUndefined).filter(noHetflag);

        for (i in atomlist) {
            if (atomlist.hasOwnProperty(i)) {
                atom = atomlist[i];

                if ((atom.atom === atomName)) {
                    if (currentChain !== atom.chain || currentResi + 1 !== atom.resi) {
                        this.drawSmoothCurve(group, points, curveWidth, colors, div);
                        points = [];
                        colors = [];
                    }
                    points.push(new TV3(atom.x, atom.y, atom.z));

                    colors.push(atom.color);
                    currentChain = atom.chain;
                    currentResi = atom.resi;
                }
            }
        }
        this.drawSmoothCurve(group, points, curveWidth, colors, div);
    };

    GLmol.prototype.drawMainchainTube = function (group, atomlist, atomName, radius) { // FIXME almost drawMainchainCurve
        var points = [],
            colors = [],
            radii = [],
            currentChain,
            currentResi,
            i,
            atom;

        atomlist = atomlist.filter(isNotUndefined).filter(noHetflag);

        for (i in atomlist) {
            if (atomlist.hasOwnProperty(i)) {
                atom = atomlist[i];

                if (atom.atom === atomName) {
                    if (currentChain !== atom.chain || currentResi + 1 !== atom.resi) {
                        this.drawSmoothTube(group, points, colors, radii);
                        points = [];
                        colors = [];
                        radii = [];
                    }
                    points.push(new TV3(atom.x, atom.y, atom.z));

                    radius = radius || (atom.b > 0) ? atom.b / 100 : 0.3;
                    radii.push(radius);
                    //if (!radius) {
                        //radii.push((atom.b > 0) ? atom.b / 100 : 0.3);
                    //} else {
                        //radii.push(radius);
                    //}

                    colors.push(atom.color);
                    currentChain = atom.chain;
                    currentResi = atom.resi;
                }
            }
        }
        this.drawSmoothTube(group, points, colors, radii);
    };

    GLmol.prototype.drawStrip = function (group, p1, p2, colors, div, thickness) {
        if ((p1.length) < 2) {
            return;
        }

        div = div || this.axisDIV;

        p1 = this.subdivide(p1, div);
        p2 = this.subdivide(p2, div);
        if (!thickness) {
            return this.drawThinStrip(group, p1, p2, colors, div);
        }

        var geo = new THREE.Geometry(),
            axis,
            p1v,
            p2v,
            a1v,
            a2v,
            i,
            j,
            lim,
            toNext,
            toSide,
            faces = [[0, 2, -6, -8], [-4, -2, 6, 4], [7, 3, -5, -1], [-3, -7, 1, 5]],
            offset,
            color,
            f,
            vsize,
            mat,
            mesh;

        for (i = 0, lim = p1.length; i < lim; i++) {
            geo.vertices.push(p1v = p1[i]); // 0
            geo.vertices.push(p1v); // 1
            geo.vertices.push(p2v = p2[i]); // 2
            geo.vertices.push(p2v); // 3
            if (i < lim - 1) {
                toNext = p1[i + 1].clone().subSelf(p1[i]);
                toSide = p2[i].clone().subSelf(p1[i]);
                axis = toSide.crossSelf(toNext).normalize().multiplyScalar(thickness);
            }
            geo.vertices.push(a1v = p1[i].clone().addSelf(axis)); // 4
            geo.vertices.push(a1v); // 5
            geo.vertices.push(a2v = p2[i].clone().addSelf(axis)); // 6
            geo.vertices.push(a2v); // 7
        }
        for (i = 1, lim = p1.length; i < lim; i++) {
            offset = 8 * i;
            color = new TCo(colors[Math.round((i - 1) / div)]);
            for (j = 0; j < 4; j++) {
                f = new THREE.Face4(offset + faces[j][0], offset + faces[j][1], offset + faces[j][2], offset + faces[j][3], undefined, color);
                geo.faces.push(f);
            }
        }
        vsize = geo.vertices.length - 8; // Cap
        for (i = 0; i < 4; i++) {
            geo.vertices.push(geo.vertices[i * 2]);
            geo.vertices.push(geo.vertices[vsize + i * 2]);
        }
        vsize += 8;
        geo.faces.push(new THREE.Face4(vsize, vsize + 2, vsize + 6, vsize + 4, undefined, geo.faces[0].color));
        geo.faces.push(new THREE.Face4(vsize + 1, vsize + 5, vsize + 7, vsize + 3, undefined, geo.faces[geo.faces.length - 3].color));

        geo.computeFaceNormals();
        geo.computeVertexNormals(false);

        mesh = this.getLambertMesh(geo);
        group.add(mesh);
    };


    GLmol.prototype.drawThinStrip = function (group, p1, p2, colors, div) {
        var geo = new THREE.Geometry(),
            i,
            lim,
            f,
            mat,
            mesh;

        for (i = 0, lim = p1.length; i < lim; i++) {
            geo.vertices.push(p1[i]); // 2i
            geo.vertices.push(p2[i]); // 2i + 1
        }
        for (i = 1, lim = p1.length; i < lim; i++) {
            f = new THREE.Face4(2 * i, 2 * i + 1, 2 * i - 1, 2 * i - 2);
            f.color = new TCo(colors[Math.round((i - 1) / div)]);
            geo.faces.push(f);
        }

        geo.computeFaceNormals();
        geo.computeVertexNormals(false);

        mesh = this.getLambertMesh(geo);
        group.add(mesh);
    };


    GLmol.prototype.IcosahedronGeometry = function () {
        this.icosahedron = this.icosahedron || new THREE.IcosahedronGeometry(1);
        return this.icosahedron;
    };


    GLmol.prototype.drawCylinder = function (group, from, to, radius, color, cap) {
        if (!from || !to) { return; }

        color = new TCo(color);
        var midpoint = new TV3().add(from, to).multiplyScalar(0.5),
            cylinderMaterial,
            cylinder,
            m;

        if (!this.cylinderGeometry) {
            this.cylinderGeometry = new THREE.CylinderGeometry(1, 1, 1, this.cylinderQuality, 1, !cap);
            this.cylinderGeometry.faceUvs = [];
            this.faceVertexUvs = [];
        }

        cylinder = this.getLambertMesh(this.cylinderGeometry, {color: color.getHex()});
        cylinder.position = midpoint;
        cylinder.lookAt(from);
        cylinder.updateMatrix();
        cylinder.matrixAutoUpdate = false;
        m = new THREE.Matrix4().makeScale(radius, radius, from.distanceTo(to));
        m.rotateX(Math.PI / 2);
        cylinder.matrix.multiplySelf(m);
        group.add(cylinder);
    };

    // FIXME: transition!
    GLmol.prototype.drawHelixAsCylinder = function (group, atomlist, radius) {
        var start = null,
            currentChain,
            currentResi,
            others = [],
            beta = [],
            i,
            atom;

        atomlist = atomlist.filter(isNotUndefined).filter(noHetflag);
        for (i in atomlist) {
            if (atomlist.hasOwnProperty(i)) {
                atom = atomlist[i];
                if ((atom.ss !== 'h' && atom.ss !== 's') || atom.ssend || atom.ssbegin) { others.push(atom.serial); }
                if (atom.ss === 's') { beta.push(atom.serial); }
                if (atom.atom !== 'CA') { continue; }

                if (atom.ss === 'h' && atom.ssend) {
                    if (start !== null) {
                        this.drawCylinder(group, new TV3(start.x, start.y, start.z), new TV3(atom.x, atom.y, atom.z), radius, atom.color, true);
                    }
                    start = null;
                }
                currentChain = atom.chain;
                currentResi = atom.resi;
                if (start === null && atom.ss === 'h' && atom.ssbegin) {
                    start = atom;
                }
            }
        }
        if (start !== null) {
            this.drawCylinder(group, new TV3(start.x, start.y, start.z), new TV3(atom.x, atom.y, atom.z), radius, atom.color);
        }
        this.drawMainchainTube(group, others, "CA", 0.3);
        this.drawStrand(group, beta, undefined, undefined, true,  0, this.helixSheetWidth, false, this.thickness * 2);
    };

    GLmol.prototype.drawCartoon = function (group, atomlist, doNotSmoothen, thickness) {
        this.drawStrand(group, atomlist, 2, undefined, true, undefined, undefined, doNotSmoothen, thickness);
    };

    GLmol.prototype.drawStrand = function (group, atomlist, num, div, fill, coilWidth, helixSheetWidth, doNotSmoothen, thickness) {
        num = num || this.strandDIV;
        div = div || this.axisDIV;
        coilWidth = coilWidth || this.coilWidth;
        doNotSmoothen = !!doNotSmoothen;
        helixSheetWidth = helixSheetWidth || this.helixSheetWidth;

        var points = [],
            colors = [],
            currentChain,
            currentResi,
            currentCA,
            prevCO = null,
            ss = null,
            ssborder = false,
            k,
            i,
            j,
            atom,
            O,
            delta,
            v;

        for (k = 0; k < num; k++) {
            points[k] = [];
        }

        atomlist = atomlist.filter(isNotUndefined).filter(noHetflag);
        for (i in atomlist) {
            if (atomlist.hasOwnProperty(i)) {
                atom = atomlist[i];

                if (atom.atom === 'O' || atom.atom === 'CA') {
                    if (atom.atom === 'CA') {
                        if (currentChain !== atom.chain || currentResi + 1 !== atom.resi) {
                            for (j = 0; !thickness && j < num; j++) {
                                this.drawSmoothCurve(group, points[j], 1, colors, div);
                            }
                            if (fill) {
                                this.drawStrip(group, points[0], points[num - 1], colors, div, thickness);
                            }

                            points = [];
                            for (k = 0; k < num; k++) { points[k] = []; }
                            colors = [];
                            prevCO = null;
                            ss = null;
                            ssborder = false;
                        }
                        currentCA = new TV3(atom.x, atom.y, atom.z);
                        currentChain = atom.chain;
                        currentResi = atom.resi;
                        ss = atom.ss;
                        ssborder = atom.ssstart || atom.ssend;
                        colors.push(atom.color);
                    } else { // O
                        O = new TV3(atom.x, atom.y, atom.z);
                        O.subSelf(currentCA);
                        O.normalize(); // can be omitted for performance
                        O.multiplyScalar((ss === 'c') ? coilWidth : helixSheetWidth);
                        if (prevCO && O.dot(prevCO) < 0) { O.negate(); }
                        prevCO = O;
                        for (j = 0; j < num; j++) {
                            delta = -1 + 2 / (num - 1) * j;
                            v = new TV3(currentCA.x + prevCO.x * delta,
                                            currentCA.y + prevCO.y * delta, currentCA.z + prevCO.z * delta);
                            if (!doNotSmoothen && ss === 's') { v.smoothen = true; }
                            points[j].push(v);
                        }
                    }
                }
            }
        }
        for (j = 0; !thickness && j < num; j++) {
            this.drawSmoothCurve(group, points[j], 1, colors, div);
        }
        if (fill) {
            this.drawStrip(group, points[0], points[num - 1], colors, div, thickness);
        }
    };

    GLmol.prototype.drawNucleicAcidLadderSub = function (geo, lineGeo, atoms, color) {
    //        color.r *= 0.9; color.g *= 0.9; color.b *= 0.9;
        var baseFaceId = geo.vertices.length,
            i,
            j,
            lim;

        function isNotUndefined(atom) {
          // assumption: glmol actually wants truthyness
            return !!atom;
        }

        if ([atoms[0], atoms[1], atoms[2], atoms[3], atoms[4], atoms[5]].every(isNotUndefined)) {
            for (i = 0; i <= 5; i++) { geo.vertices.push(atoms[i]); }
            geo.faces.push(new TF3(baseFaceId, baseFaceId + 1, baseFaceId + 2));
            geo.faces.push(new TF3(baseFaceId, baseFaceId + 2, baseFaceId + 3));
            geo.faces.push(new TF3(baseFaceId, baseFaceId + 3, baseFaceId + 4));
            geo.faces.push(new TF3(baseFaceId, baseFaceId + 4, baseFaceId + 5));
            for (j = geo.faces.length - 4, lim = geo.faces.length; j < lim; j++) {
                geo.faces[j].color = color;
            }
        }
        if ([atoms[3], atoms[4], atoms[6], atoms[7], atoms[7], atoms[8]].every(isNotUndefined)) {
            geo.vertices.push(atoms[4]);
            geo.vertices.push(atoms[3]);
            geo.vertices.push(atoms[6]);
            geo.vertices.push(atoms[7]);
            geo.vertices.push(atoms[8]);
            for (i = 0; i <= 4; i++) {
                geo.colors.push(color);
            }
            geo.faces.push(new TF3(baseFaceId, baseFaceId + 1, baseFaceId + 2));
            geo.faces.push(new TF3(baseFaceId, baseFaceId + 2, baseFaceId + 3));
            geo.faces.push(new TF3(baseFaceId, baseFaceId + 3, baseFaceId + 4));
            for (j = geo.faces.length - 3, lim = geo.faces.length; j < lim; j++) {
                geo.faces[j].color = color;
            }
        }
    };

    GLmol.prototype.drawNucleicAcidLadder = function (group, atomlist) {
        var geo = new THREE.Geometry(),
            lineGeo = new THREE.Geometry(),
            baseAtoms = ["N1", "C2", "N3", "C4", "C5", "C6", "N9", "C8", "N7"],
            currentChain,
            currentResi,
            currentComponent = new Array(baseAtoms.length),
            color = new TCo(0xcc0000), // MAGIC
            i,
            atom,
            pos,
            mat,
            mesh;

        atomlist = atomlist.filter(isNotUndefined).filter(noHetflag);

        for (i in atomlist) {
            if (atomlist.hasOwnProperty(i)) {
                atom = atomlist[i];

                if (atom.resi !== currentResi || atom.chain !== currentChain) {
                    this.drawNucleicAcidLadderSub(geo, lineGeo, currentComponent, color);
                    currentComponent = new Array(baseAtoms.length);
                }
                pos = baseAtoms.indexOf(atom.atom);
                if (pos !== -1) {
                    currentComponent[pos] = new TV3(atom.x, atom.y, atom.z);
                }
                if (atom.atom === 'O3\'') {
                    color = new TCo(atom.color);
                }
                currentResi = atom.resi;
                currentChain = atom.chain;
            }
        }

        this.drawNucleicAcidLadderSub(geo, lineGeo, currentComponent, color);

        geo.computeFaceNormals();

        mesh = this.getLambertMesh(geo);
        group.add(mesh);
    };

    GLmol.prototype.drawNucleicAcidStick = function (group, atomlist) {
        var currentChain,
            currentResi,
            start = null,
            end = null,
            i,
            atom;

        atomlist = atomlist.filter(isNotUndefined).filter(noHetflag);

        for (i in atomlist) {
            if (atomlist.hasOwnProperty(i)) {
                atom = atomlist[i];

                if (atom.resi !== currentResi || atom.chain !== currentChain) {
                    if (start !== null && end !== null) {
                        this.drawCylinder(group, new TV3(start.x, start.y, start.z), new TV3(end.x, end.y, end.z), 0.3, start.color, true);
                    }
                    start = null;
                    end = null;
                }
                if (atom.atom === 'O3\'') {
                    start = atom;
                }
                if (atom.resn === '  A' || atom.resn === '  G' || atom.resn === ' DA' || atom.resn === ' DG') {
                    if (atom.atom === 'N1') {
                        end = atom; //  N1(AG), N3(CTU)
                    }
                } else if (atom.atom === 'N3') {
                    end = atom;
                }
                currentResi = atom.resi;
                currentChain = atom.chain;
            }
        }
        if (start !== null && end !== null) {
            this.drawCylinder(group, new TV3(start.x, start.y, start.z), new TV3(end.x, end.y, end.z), 0.3, start.color, true);
        }
    };

    GLmol.prototype.drawNucleicAcidLine = function (group, atomlist) {
        var currentChain,
            currentResi,
            start = null,
            end = null,
            geo = new THREE.Geometry(),
            i,
            atom,
            mat,
            line;

        atomlist = atomlist.filter(isNotUndefined).filter(noHetflag);

        for (i in atomlist) {
            if (atomlist.hasOwnProperty(i)) {

                atom = atomlist[i];

                if (atom.resi !== currentResi || atom.chain !== currentChain) {
                    if (start !== null && end !== null) {
                        geo.vertices.push(new TV3(start.x, start.y, start.z));
                        geo.colors.push(new TCo(start.color));
                        geo.vertices.push(new TV3(end.x, end.y, end.z));
                        geo.colors.push(new TCo(start.color));
                    }
                    start = null;
                    end = null;
                }
                if (atom.atom === 'O3\'') {
                    start = atom;
                }
                if (atom.resn === '  A' || atom.resn === '  G' || atom.resn === ' DA' || atom.resn === ' DG') {
                    if (atom.atom === 'N1') {
                        end = atom; //  N1(AG), N3(CTU)
                    }
                } else if (atom.atom === 'N3') {
                    end = atom;
                }
                currentResi = atom.resi;
                currentChain = atom.chain;
            }
        }
        if (start !== null && end !== null) {
            geo.vertices.push(new TV3(start.x, start.y, start.z));
            geo.colors.push(new TCo(start.color));
            geo.vertices.push(new TV3(end.x, end.y, end.z));
            geo.colors.push(new TCo(start.color));
        }

        mat =  new THREE.LineBasicMaterial({linewidth: 1, linejoin: false});
        mat.linewidth = 1.5;
        mat.vertexColors = true;
        line = new THREE.Line(geo, mat, THREE.LinePieces);
        group.add(line);
    };

    GLmol.prototype.drawCartoonNucleicAcid = function (group, atomlist, div, thickness) {
        this.drawStrandNucleicAcid(group, atomlist, 2, div, true, undefined, thickness);
    };

    GLmol.prototype.drawStrandNucleicAcid = function (group, atomlist, num, div, fill, nucleicAcidWidth, thickness) {
        nucleicAcidWidth = nucleicAcidWidth || this.nucleicAcidWidth;
        div = div || this.axisDIV;
        num = num || this.nucleicAcidStrandDIV;

        var points = [],
            colors = [],
            currentChain,
            currentResi,
            currentO3,
            prevOO = null,
            i,
            k,
            j,
            atom,
            O,
            delta;

        for (k = 0; k < num; k++) {
            points[k] = [];
        }

        atomlist = atomlist.filter(isNotUndefined).filter(noHetflag);

        for (i in atomlist) {
            if (atomlist.hasOwnProperty(i)) {
                atom = atomlist[i];

                if (atom.atom === 'O3\'' || atom.atom === 'OP2') {
                    if (atom.atom === 'O3\'') { // to connect 3' end. FIXME: better way to do?
                        if (currentChain !== atom.chain || currentResi + 1 !== atom.resi) {
                            if (currentO3) {
                                for (j = 0; j < num; j++) {
                                    delta = -1 + 2 / (num - 1) * j;
                                    points[j].push(new TV3(currentO3.x + prevOO.x * delta, currentO3.y + prevOO.y * delta, currentO3.z + prevOO.z * delta));
                                }
                            }
                            if (fill) {
                                this.drawStrip(group, points[0], points[1], colors, div, thickness);
                            }
                            for (j = 0; !thickness && j < num; j++) {
                                this.drawSmoothCurve(group, points[j], 1, colors, div);
                            }

                            points = [];
                            for (k = 0; k < num; k++) {
                                points[k] = [];
                            }
                            colors = [];
                            prevOO = null;
                        }
                        currentO3 = new TV3(atom.x, atom.y, atom.z);
                        currentChain = atom.chain;
                        currentResi = atom.resi;
                        colors.push(atom.color);
                    } else { // OP2
                        if (!currentO3) { // for 5' phosphate (e.g. 3QX3)
                            prevOO = null;
                            continue;
                        }
                        O = new TV3(atom.x, atom.y, atom.z);
                        O.subSelf(currentO3);
                        O.normalize().multiplyScalar(nucleicAcidWidth);  // TODO: refactor
                        if (prevOO && O.dot(prevOO) < 0) {
                            O.negate();
                        }
                        prevOO = O;
                        for (j = 0; j < num; j++) {
                            delta = -1 + 2 / (num - 1) * j;
                            points[j].push(new TV3(currentO3.x + prevOO.x * delta, currentO3.y + prevOO.y * delta, currentO3.z + prevOO.z * delta));
                        }
                        currentO3 = null;
                    }
                }
            }
        }
        if (currentO3) {
            for (j = 0; j < num; j++) {
                delta = -1 + 2 / (num - 1) * j;
                points[j].push(new TV3(currentO3.x + prevOO.x * delta, currentO3.y + prevOO.y * delta, currentO3.z + prevOO.z * delta));
            }
        }
        if (fill) {
            this.drawStrip(group, points[0], points[1], colors, div, thickness);
        }
        for (j = 0; !thickness && j < num; j++) {
            this.drawSmoothCurve(group, points[j], 1, colors, div);
        }
    };

    GLmol.prototype.drawDottedLines = function (group, points, color) {
        var geo = new THREE.Geometry(),
            step = 0.3,
            i,
            lim,
            delta,
            dist,
            jlim,
            j,
            p,
            p1,
            p2,
            mat,
            line;

        for (i = 0, lim = Math.floor(points.length / 2); i < lim; i++) {
            p1 = points[2 * i];
            p2 = points[2 * i + 1];
            delta = p2.clone().subSelf(p1);
            dist = delta.length();
            delta.normalize().multiplyScalar(step);
            jlim =  Math.floor(dist / step);
            for (j = 0; j < jlim; j++) {
                p = new TV3(p1.x + delta.x * j, p1.y + delta.y * j, p1.z + delta.z * j);
                geo.vertices.push(p);
            }
            if (jlim % 2 === 1) {
                geo.vertices.push(p2);
            }
        }

        mat = new THREE.LineBasicMaterial({'color': color.getHex()});
        mat.linewidth = 2;
        line = new THREE.Line(geo, mat, THREE.LinePieces);
        group.add(line);
    };



    /* helper functions, generalize them! */
    function isNotUndefined(atom) {
        return !!atom;
    }
    function hasHetflag(atom) {
        return !!atom.hetflag;
    }

    function noHetflag(atom) {
        return hasHetflag(atom) === false;
    }

    function notCA(atom) {
        return atom.atom !== 'CA';
    }

    GLmol.isNotSolvent = function (atom) {
        return atom.resn !== 'HOH';
    };

    isNotSolvent = GLmol.isNotSolvent;

    function getAtomSerial(atom) {
        console.warn("getAtomSerial is deprecated, use atom's index instead");
        return atom;
        //return atom.serial;
    }

    /* */

    GLmol.prototype.getAllAtoms = function () {
        console.warn("getAllAtoms is deprecated, use .atoms instead");
        return this.atoms.map(getAtomSerial);
    };

    GLmol.prototype.getAtoms = function (atomlist) {
        console.warn("getAtoms is deprecated, access atoms list directly");
        return atomlist;
        //return this.atoms.filter(function (atom, index) {
            //if (atomlist.indexOf(index) !== -1) { return true; }
        //});
    };


    GLmol.prototype.getHetatms = function (atomlist) {
        return atomlist.filter(isNotUndefined).filter(hasHetflag);
    };

    GLmol.prototype.removeSolvents = function (atomlist) {
        return atomlist.filter(isNotUndefined).filter(isNotSolvent);
    };

    GLmol.prototype.getProteins = function (atomlist) {
        return atomlist.filter(isNotUndefined).filter(noHetflag);
    };

    // TODO: Test
    GLmol.prototype.excludeAtoms = function (atomlist, deleteList) {
        return atomlist.filter(function (atom) { return (deleteList.indexOf(atom) === -1); });
    };


    GLmol.prototype.getSidechains = function (atomlist) {
        //TODO: what is this called?
        function isSidechain(atom) {
            return !(atom.atom === 'C' || atom.atom === 'O' || (atom.atom === 'N' && atom.resn !== "PRO"));
        }
        return atomlist.filter(isNotUndefined).filter(noHetflag).filter(isSidechain);
    };

    GLmol.prototype.getAtomsWithin = function (atomlist, extent) {
        return atomlist.filter(isNotUndefined).filter(function (atom) {
            if (atom.x < extent[0][0] || atom.x > extent[1][0]) { return false; }
            if (atom.y < extent[0][1] || atom.y > extent[1][1]) { return false; }
            if (atom.z < extent[0][2] || atom.z > extent[1][2]) { return false; }

            return true;
        });
    };

    GLmol.prototype.getExtent = function (atomlist) {
        var xmin,
            ymin,
            zmin,
            xmax,
            ymax,
            zmax,
            xsum,
            ysum,
            zsum,
            cnt,
            atom,
            i;

        xmin = ymin = zmin = 9999;
        xmax = ymax = zmax = -999;
        xsum = ysum = zsum = cnt = 0;

        atomlist = atomlist.filter(isNotUndefined);
        for (i in atomlist) {
            if (atomlist.hasOwnProperty(i)) {
                atom = atomlist[i];
                cnt++;

                xsum += atom.x;
                ysum += atom.y;
                zsum += atom.z;

                xmin = (xmin < atom.x) ? xmin : atom.x;
                ymin = (ymin < atom.y) ? ymin : atom.y;
                zmin = (zmin < atom.z) ? zmin : atom.z;
                xmax = (xmax > atom.x) ? xmax : atom.x;
                ymax = (ymax > atom.y) ? ymax : atom.y;
                zmax = (zmax > atom.z) ? zmax : atom.z;
            }
        }
        return [[xmin, ymin, zmin], [xmax, ymax, zmax], [xsum / cnt, ysum / cnt, zsum / cnt]];
    };

    GLmol.prototype.getResiduesById = function (atomlist, resi) { // FIXME: this is almost the same as excludeAtoms
        return atomlist.filter(function (atom) { return (resi.indexOf(atom.resi) === -1); });
    };

    GLmol.prototype.getResidueBySS = function (atomlist, ss) { // FIXME and again.
        return atomlist.filter(function (atom) { return (ss.indexOf(atom.ss) === -1); });
    };

    GLmol.prototype.getChain = function (atomlist, chain) {
        var chains = {},
            i,
            lim,
            atom;

        console.log(chain);
        chain = chain.toString(); // concat if Array

        for (i = 0, lim = chain.length; i < lim; i++) {
            chains[chain.substr(i, 1)] = true;
        }

        return atomlist.filter(function (atom) {
            return chains[atom.chain];
        });
    };

    // for HETATM only
    GLmol.prototype.getNonbonded = function (atomlist, chain) { // XXX chain argument unused!!
        return atomlist.filter(isNotUndefined).filter(hasHetflag).filter(
            function (atom) { return atom.bonds.length === 0; }
        );
    };

    GLmol.prototype.colorByAtom = function (atomlist, colors) {
        atomlist.filter(isNotUndefined).forEach(function (atom) {
            atom.color = this.ElementColors[atom.elem] || this.defaultColor || colors[atom.elem];
        }.bind(this));
    };


    // MEMO: Color only CA. maybe I should add atom.cartoonColor.
    GLmol.prototype.colorByStructure = function (atomlist, helixColor, sheetColor, colorSidechains) {
        atomlist.filter(isNotUndefined).filter(noHetflag).filter(notCA).forEach(function (atom) {
            if (atom.ss[0] === 's') {
                atom.color = sheetColor;
            } else if (atom.ss[0] === 'h') {
                atom.color = helixColor;
            }
        });
    };


    GLmol.prototype.colorByBFactor = function (atomlist, colorSidechains) {
        var minB = 1000, // XXX magic values
            maxB = -1000,
            i,
            atom,
            mid = (maxB + minB) / 2,
            range = (maxB - minB) / 2,
            color;

        atomlist = atomlist.filter(isNotUndefined).filter(noHetflag);
        for (i in atomlist) {
            if (atomlist.hasOwnProperty(i)) {
                atom = atomlist[i];
                if (colorSidechains || atom.atom === 'CA' || atom.atom === 'O3\'') {
                    if (minB > atom.b) { minB = atom.b; }
                    if (maxB < atom.b) { maxB = atom.b; }
                }
            }
        }

        if (range < 0.01 && range > -0.01) { return; }
        for (i in atomlist) {
            if (atomlist.hasOwnProperty(i)) {
                atom = atomlist[i];
                if (colorSidechains || atom.atom === 'CA' || atom.atom === 'O3\'') {
                    color = new TCo(0);
                    if (atom.b < mid) {
                        color.setHSV(0.667, (mid - atom.b) / range, 1);
                    } else {
                        color.setHSV(0, (atom.b - mid) / range, 1);
                        atom.color = color.getHex();
                    }
                }
            }
        }
    };

    GLmol.prototype.colorByChain = function (atomlist, colorSidechains) {
        var i, atom, color;

        atomlist = atomlist.filter(isNotUndefined).filter(noHetflag);
        for (i in atomlist) {
            if (atomlist.hasOwnProperty(i)) {
                atom = atomlist[i];
                if (colorSidechains || atom.atom === 'CA' || atom.atom === 'O3\'') {
                    color = new TCo(0);
                    color.setHSV((atom.chain.charCodeAt(0) * 5) % 17 / 17.0, 1, 0.9);
                    atom.color = color.getHex();
                }
            }
        }
    };

    GLmol.prototype.colorByResidue = function (atomlist, residueColors) {
        atomlist.filter(isNotUndefined).forEach(
            function (atom) {
                var c = residueColors[atom.resn];
                if (c) {
                    atom.color = c;
                }
            }
        );
    };

    GLmol.prototype.colorAtoms = function (atomlist, c) {
        atomlist.filter(isNotUndefined).forEach(
            function (atom) {
                atom.color = c;
            }
        );
    };

    GLmol.prototype.colorByPolarity = function (atomlist, polar, nonpolar) {
        var polarResidues = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS'],
            nonPolarResidues = ['GLY', 'PRO', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TYR', 'TRP'],
            colorMap = {},
            i;

        for (i in polarResidues) {
            if (polarResidues.hasOwnProperty(i)) {
                colorMap[polarResidues[i]] = polar;
            }
        }
        for (i in nonPolarResidues) {
            if (nonPolarResidues.hasOwnProperty(i)) {
                colorMap[nonPolarResidues[i]] = nonpolar;
            }
        }
        this.colorByResidue(atomlist, colorMap);
    };

    // TODO: Add near(atomlist, neighbor, distanceCutoff)
    // TODO: Add expandToResidue(atomlist)

    GLmol.prototype.colorChainbow = function (atomlist, colorSidechains) {
        var cnt = 0,
            atom,
            i,
            total,
            color;

        atomlist = atomlist.filter(isNotUndefined).filter(noHetflag);

        for (i in atomlist) {
            if (atomlist.hasOwnProperty(i)) {
                atom = atomlist[i];

                if (colorSidechains || atom.atom !== 'CA' || atom.atom !== 'O3\'') {
                    cnt++;
                }
            }
        }

        total = cnt;
        cnt = 0;

        for (i in atomlist) {
            if (atomlist.hasOwnProperty(i)) {
                atom = atomlist[i];

                if (colorSidechains || atom.atom !== 'CA' || atom.atom !== 'O3\'') {
                    color = new TCo(0);
                    color.setHSV(240.0 / 360 * (1 - cnt / total), 1, 0.9);
                    atom.color = color.getHex();
                    cnt++;
                }
            }
        }
    };

    GLmol.prototype.drawSymmetryMates2 = function (group, asu, matrices) {
        if (!matrices) { return; }
        asu.matrixAutoUpdate = false;

        var cnt = 1,
            i,
            j,
            mat,
            symmetryMate;

        this.protein.appliedMatrix = new THREE.Matrix4();
        for (i = 0; i < matrices.length; i++) {
            mat = matrices[i];
            if (!mat || mat.isIdentity()) { continue; }
            //console.log(mat);
            symmetryMate = THREE.SceneUtils.cloneObject(asu);
            symmetryMate.matrix = mat;
            group.add(symmetryMate);
            for (j = 0; j < 16; j++) {
                this.protein.appliedMatrix.elements[j] += mat.elements[j];
            }
            cnt++;
        }
        this.protein.appliedMatrix.multiplyScalar(cnt);
    };


    GLmol.prototype.drawSymmetryMatesWithTranslation2 = function (group, asu, matrices) {
        if (!matrices) { return; }
        var p = this.protein,
            i,
            a,
            b,
            c,
            mat,
            translationMat,
            symop,
            symmetryMate;
        asu.matrixAutoUpdate = false;

        for (i = 0; i < matrices.length; i++) {
            mat = matrices[i];
            if (!mat) { continue; }

            for (a = -1; a <= 0; a++) {
                for (b = -1; b <= 0; b++) {
                    for (c = -1; c <= 0; c++) {
                        translationMat = new THREE.Matrix4().makeTranslation(
                            p.ax * a + p.bx * b + p.cx * c,
                            p.ay * a + p.by * b + p.cy * c,
                            p.az * a + p.bz * b + p.cz * c
                        );
                        symop = mat.clone().multiplySelf(translationMat);
                        if (symop.isIdentity()) { continue; }
                        symmetryMate = THREE.SceneUtils.cloneObject(asu);
                        symmetryMate.matrix = symop;
                        group.add(symmetryMate);
                    }
                }
            }
        }
    };

    GLmol.prototype.createTextTex = function (text, size, color) {
        var canvas = document.createElement("canvas"),
            ctx = canvas.getContext("2d"),
            tex;

        canvas.style.backgroundColor = "rgba(0, 0, 0, 0.0)";
        ctx.font = size + "pt Arial";
        canvas.width = ctx.measureText(text).width;
        canvas.height = size; //size; // This resets fonts, so we have to set it again
        ctx.fillStyle = color || "rgba(0, 0, 0, 1.0)";
        ctx.strokeStyle = ctx.fillStyle;
        ctx.font = size + "pt Arial";
        ctx.fillText(text, 0, size * 0.9);

        tex = new THREE.Texture(canvas);
        tex.needsUpdate = true;
        tex.magFilter = tex.minFilter = THREE.LinearFilter;
        return tex;
    };

    GLmol.prototype.billboard = function (texture) {
        var sprite = new THREE.Sprite({
            map: texture,
            useScreenCoordinates: false,
            transparent: true
        });
        return sprite;
    };

    GLmol.prototype.labelAtom = function (atom, text, show) {
        var texture = this.createTextTex(text, 40, "#fff");
        var bb = this.billboard(texture);
        bb.scale.set( texture.image.width / 40 , texture.image.height / 40, 1.0 )



        bb.position.set(atom.x, atom.y, atom.z);
        this.modelGroup.add(bb);

        atom.label = bb;

        if (show) {
            this.show();
        }
    };

    GLmol.prototype.defineRepresentation = function (all) {
        var all = all || this.atoms,
            hetatm = this.getHetatms(all).filter(isNotSolvent);

        var options = this.options || {
            colorBy: "chain",
            mainChainAs: "thick ribbon",
            doNotSmoothen: false,
            sideChainAsLines: "true",
            nonBondedAs: "stars",
            hetatmsAs: "ball and stick",
            // --
            //
            //unitCell: true,
            //biologicalAssembly: true,
            //crystalPacking: true,
            //labelCA: true,
        };

        var asu = new THREE.Object3D();

        this.colorByAtom(all, {});

        switch (options.colorBy) {
        case "spectrum":
            this.colorChainbow(all);
            break;
        case "secondary structure":
            this.colorByStructure(all, 0xcc00cc, 0x00cccc);
            break;
        case "B factor":
            this.colorByBFactor(all);
            break;
        case "polar/nonpolar":
            this.colorByPolarity(all, 0xcc0000, 0xcccccc);
            break;
        case "chain":
        default:
            this.colorByChain(all);
            break;
        }

        var doNotSmoothen = options.doNotSmoothen || false;

        switch (options.mainChainAs) {
        case "thin ribbon and lines":
            this.drawCartoon(asu, all, doNotSmoothen);
            this.drawCartoonNucleicAcid(asu, all);
            this.drawBondsAsLine(asu, all, this.lineWidth);
            break;
        case "thin ribbon":
            this.drawCartoon(asu, all, doNotSmoothen);
            this.drawCartoonNucleicAcid(asu, all);
            break;
        case "strand":
            this.drawStrand(asu, all, null, null, null, null, null, doNotSmoothen);
            this.drawStrandNucleicAcid(asu, all);
            break;
        case "cylinder and plate":
            this.drawHelixAsCylinder(asu, all, 1.6);
            this.drawCartoonNucleicAcid(asu, all);
            break;
        case "C alpha trace":
            this.drawMainchainCurve(asu, all, this.curveWidth, 'CA', 1);
            this.drawMainchainCurve(asu, all, this.curveWidth, 'O3\'', 1);
            break;
        case "B factor Tube":
            this.drawMainchainTube(asu, all, 'CA');
            this.drawMainchainTube(asu, all, 'O3\''); // FIXME: 5' end problem!
            break;
        case "bonds":
            this.drawBondsAsLine(asu, all, this.lineWidth);
            break;
        case "thick ribbon":
        default:
            this.drawCartoon(asu, all, doNotSmoothen, this.thickness);
            this.drawCartoonNucleicAcid(asu, all, null, this.thickness);
            break;
        }

        if (options.sideChainAsLines) {
            this.drawBondsAsLine(asu, this.getSidechains(all), this.lineWidth);
        }

        
        if (options.mainChainAs) {
            var allHet = this.getHetatms(all);
            var nonBonded = this.getNonbonded(allHet);
            switch (options.mainChainAs) {
            case "stars":
                this.drawAsCross(target, nonBonded, 0.3);
                break;
            case "spheres":
                this.drawAtomsAsIcosahedron(target, nonBonded, 0.3, true);
                break;
            }


        }

        switch (options.hetatmsAs) {
            case "sticks":
             this.drawBondsAsStick(asu, hetatm, this.cylinderRadius, this.cylinderRadius, true); //hetatm sticks
             this.drawCartoon(asu, all, this.curveWidth, this.thickness);
                break;
            case "ball and stick":
                this.drawBondsAsStick(asu, hetatm, this.cylinderRadius / 2.0, this.cylinderRadius, true, false, 0.3);
                break;
            case "ball and stick multiple":
                this.drawBondsAsStick(asu, hetatm, this.cylinderRadius / 2.0, this.cylinderRadius, true, true, 0.3);
                break;
            case "spheres":
                this.drawAtomsAsSphere(asu, hetatm, this.sphereRadius);
                break;
            case "icosahedrons":
                this.drawAtomsAsIcosahedron(asu, hetatm, this.sphereRadius);
                break;
            case "lines":
                this.drawBondsAsLine(asu, hetatm, this.curveWidth);
                break;
        }

        if (options.unitCell) {
            this.drawUnitcell(this.modelGroup);
        }

        if (options.biologicalAssembly) {
            this.drawSymmetryMates2(this.modelGroup, asu, this.protein.biomtMatrices);
        }

        if (options.crystalPacking) {
            this.drawSymmetryMatesWithTranslation2(this.modelGroup, asu, this.protein.symMat);
        }
        if (options.labelCA) {
            console.time("labelCA");
            this.atoms.filter(function(atom){return atom.atom == "CA"}).forEach(function(CA) {
                this.labelAtom(CA, CA.resn + "" + CA.resi, false);
            }.bind(this))
            this.show();
            console.timeEnd("labelCA");
        }

        this.modelGroup.add(asu);

        if (options.proteinSurface) {
            console.log("generating proteinsurface");
            this.generateMesh(this.modelGroup, this.atoms, options.proteinSurface);
        }


    };

    GLmol.prototype.getView = function () {
        if (!this.modelGroup) {
            return [0, 0, 0, 0, 0, 0, 0, 1];
        }

        var pos = this.modelGroup.position,
            q = this.rotationGroup.quaternion;
        return [pos.x, pos.y, pos.z, this.rotationGroup.position.z, q.x, q.y, q.z, q.w];
    };

    GLmol.prototype.setView = function (arg) {
        if (!this.modelGroup || !this.rotationGroup) { return; }
        this.modelGroup.position.x = arg[0];
        this.modelGroup.position.y = arg[1];
        this.modelGroup.position.z = arg[2];
        this.rotationGroup.position.z = arg[3];
        this.rotationGroup.quaternion.x = arg[4];
        this.rotationGroup.quaternion.y = arg[5];
        this.rotationGroup.quaternion.z = arg[6];
        this.rotationGroup.quaternion.w = arg[7];
        this.show();
    };

    GLmol.prototype.setBackground = function (hex, a) {
        a = a || 1.0;
        this.bgColor = hex;
        this.renderer.setClearColorHex(hex, a);
        this.scene.fog.color = new TCo(hex);
    };

    GLmol.prototype.initializeScene = function () {
        // CHECK: Should I explicitly call scene.deallocateObject?
        this.scene = new THREE.Scene();
        this.scene.fog = new THREE.Fog(this.bgColor, 100, 200);

        this.modelGroup = new THREE.Object3D();
        this.rotationGroup = new THREE.Object3D();
        this.rotationGroup.useQuaternion = true;
        this.rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0);
        this.rotationGroup.add(this.modelGroup);

        this.scene.add(this.rotationGroup);
        this.setupLights(this.scene);
    };

    GLmol.prototype.rebuildScene = function () {
        console.time("built scene");
        this.initializeScene();
        this.defineRepresentation();
        //var view = this.getView(); // seems unnecessary to get and then set the same view.
        //this.setView(view);
        this.show();
        console.timeEnd("built scene")

    };

    GLmol.prototype.zoomInto = function (atomlist, keepSlab) {
        var tmp = this.getExtent(atomlist),
            center = new TV3(tmp[2][0], tmp[2][1], tmp[2][2]),
            x,
            y,
            z,
            maxD;

        if (this.protein.appliedMatrix) {
            center = this.protein.appliedMatrix.multiplyVector3(center);
        }

        this.modelGroup.position = center.multiplyScalar(-1);

        x = tmp[1][0] - tmp[0][0];
        y = tmp[1][1] - tmp[0][1];
        z = tmp[1][2] - tmp[0][2];

        maxD = Math.sqrt(x * x + y * y + z * z);
        maxD = Math.max(maxD, 25);

        if (!keepSlab) {
            this.slabNear = -maxD / 1.9;
            this.slabFar = maxD / 3;
        }

        this.rotationGroup.position.z = maxD * 0.35 / Math.tan(Math.PI / 180.0 * this.camera.fov / 2) - 150;
        this.rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0);
    };


    GLmol.prototype.loadMolecule = function (repressZoom) {
        //var source = document.querySelector(this.queryselector + '_src').innerHTML;
        var source = $(this.queryselector + '_src').val();
        this.loadMoleculeStr(repressZoom, source);
    };

    GLmol.prototype.loadMoleculeStr = function (repressZoom, source) {
        console.time("parsed")
        var title_elem = document.querySelector(this.queryselector + '_pdbTitle'),
            titleStr = '',
            parser_id,
            parsers = [this.parseSDF, this.parseXYZ, this.parsePDB2];


        for (parser_id in parsers) {
            if (parsers.hasOwnProperty(parser_id)) {
                if (parsers[parser_id].call(this, source)) {
                    break;
                }
            }
        }
        if (this.atoms.length === 0) {
          throw Error("No atoms found in parsed file, must be an error");
        }
        console.timeEnd("parsed")

        if (this.protein.pdbID) {
            titleStr += '<a href="http://www.rcsb.org/pdb/explore/explore.do?structureId=' + this.protein.pdbID + '">' + this.protein.pdbID + '</a>';
        }
        if (this.protein.title) {
            titleStr += '<br>' + this.protein.title;
        }

        if (title_elem) {
            title_elem.innerHTML = titleStr;
        }// jQ's method is more thorough

        this.rebuildScene();

        if (!repressZoom) {
            this.zoomInto(this.atoms);
        }

        this.show();
    };

    GLmol.prototype.setSlabAndFog = function () {
        var center = this.rotationGroup.position.z - this.camera.position.z;
        if (center < 1) {
            center = 1;
        }
        this.camera.near = center + this.slabNear;
        this.camera.near = Math.max(this.camera.near, 1);

        this.camera.far = center + this.slabFar;
        if (this.camera.near + 1 > this.camera.far) {
            this.camera.far = this.camera.near + 1;
        }

        if (this.camera instanceof THREE.PerspectiveCamera) {
            this.camera.fov = this.fov;
        } else {
            this.camera.right = center * Math.tan(Math.PI / 180 * this.fov);
            this.camera.left = -this.camera.right;
            this.camera.top = this.camera.right / this.ASPECT;
            this.camera.bottom = -this.camera.top;
        }
        this.camera.updateProjectionMatrix();
        this.scene.fog.near = this.camera.near + this.fogStart * (this.camera.far - this.camera.near);
        //   if (this.scene.fog.near > center) this.scene.fog.near = center;
        this.scene.fog.far = this.camera.far;
    };

    GLmol.prototype.adjustPos = function (ev) {
        var x = ev.pageX, y = ev.pageY;
        if (!ev.originalEvent) {
          // Zepto compatibility
          ev.originalEvent = ev;
        }
        
        if (ev.originalEvent.targetTouches && ev.originalEvent.targetTouches[0]) {
            x = ev.originalEvent.targetTouches[0].pageX;
            y = ev.originalEvent.targetTouches[0].pageY;
        }
        ev.x = x;
        ev.y = y;
    };

    GLmol.prototype.enableMouse = function () {
        var glDOM = $(this.renderer.domElement);

       // TODO: Better touch panel support.
       // Contribution is needed as I don't own any iOS or Android device with WebGL support.
        glDOM.on('mousedown touchstart', function (ev) {
            ev.preventDefault();
            if (!this.scene) { return; }
            this.adjustPos(ev);
            if (!ev.x) { return; }
            this.isDragging = true;
            this.mouseButton = ev.which;
            this.mouseStartX = ev.x;
            this.mouseStartY = ev.y;
            this.cq = this.rotationGroup.quaternion;
            this.cz = this.rotationGroup.position.z;
            this.currentModelPos = this.modelGroup.position.clone();
            this.cslabNear = this.slabNear;
            this.cslabFar = this.slabFar;
            glDOM.parent().addClass("active")
        }.bind(this));

        glDOM.on('DOMMouseScroll mousewheel', function (ev) { // Zoom
            ev.preventDefault();
            if (!ev.originalEvent) {
              // Zepto compatibility
              ev.originalEvent = ev;
            }
            if (!this.scene) { return; }
            var scaleFactor = (this.rotationGroup.position.z - this.CAMERA_Z) * 0.85;
            if (ev.originalEvent.detail) { // Webkit
                this.rotationGroup.position.z += scaleFactor * ev.originalEvent.detail / 10;
            } else if (ev.originalEvent.wheelDelta) { // Firefox
                this.rotationGroup.position.z -= scaleFactor * ev.originalEvent.wheelDelta / 400;
            }
            console.log(ev.originalEvent.wheelDelta, ev.originalEvent.detail, this.rotationGroup.position.z);
            this.show();
        }.bind(this));

        glDOM.on("contextmenu", function (ev) { ev.preventDefault(); });

        $("body").on('mouseup touchend', function (ev) {
            glDOM.parent().removeClass("active")
            this.isDragging = false;
            var x,
                y,
                dx,
                dy,
                r,
                mvMat,
                pmvMat,
                pmvMatInv,
                nearest,
                i,
                atom,
                v,
                r2,
                tx,
                ty,
                ilim,
                bb;
            this.isDragging = false;

            this.adjustPos(ev);
            x = ev.x;
            y = ev.y;
            if (!x) {
                return;
            }
            dx = x - this.mouseStartX;
            dy = y - this.mouseStartY;
            r = Math.sqrt(dx * dx + dy * dy);
            if (r > 2) {
                return;
            }
            x -= this.container.offset().left //Zepto fix, originally this.container.position().left;
            y -= this.container.offset().top //Zepto fix, originally this.container.position().top;


            mvMat = new THREE.Matrix4().multiply(this.camera.matrixWorldInverse, this.modelGroup.matrixWorld);
            pmvMat = new THREE.Matrix4().multiply(this.camera.projectionMatrix, mvMat);
            pmvMatInv = new THREE.Matrix4().getInverse(pmvMat);
            tx = x / this.WIDTH * 2 - 1;
            ty = 1 - y / this.HEIGHT * 2;
            nearest = [1, {}, new TV3(0, 0, 1000)];
            for (i = 0, ilim = this.atoms.length; i < ilim; i++) {
                atom = this.atoms[i];
                if (!atom) { continue; }
                if (atom.resn === "HOH") { continue; }

                v = new TV3(atom.x, atom.y, atom.z);
                pmvMat.multiplyVector3(v);
                r2 = (v.x - tx) * (v.x - tx) + (v.y - ty) * (v.y - ty);
                if (r2 > 0.0005) { continue; }
                if (v.z < nearest[2].z) { nearest = [r2, atom, v]; }
                if (r2 > 0.0002) { continue; }
                if (r2 < nearest[0]) { nearest = [r2, atom, v]; }
            }
            atom = nearest[1];
            if (Object.keys(atom).length === 0) { return; }

            glDOM.trigger("atom:select", [atom]);

        }.bind(this));

        glDOM.on('mousemove touchmove', function (ev) { // touchmove
            var mode = this.mouseMode || 0,
                //modeRadio = document.querySelectorAll('input[name=' + this.id + '_mouseMode]:checked'),
                dx,
                dy,
                r,
                scaleFactor,
                translationByScreen,
                q,
                qinv,
                translation,
                rs,
                x,
                y;

            ev.preventDefault();
            if (!this.scene) { return; }
            if (!this.isDragging) { return; }

            this.adjustPos(ev);
            x = ev.x;
            y = ev.y;
            if (!x) { return; }
            dx = (x - this.mouseStartX) / this.WIDTH;
            dy = (y - this.mouseStartY) / this.HEIGHT;
            r = Math.sqrt(dx * dx + dy * dy);
            if (mode === 3 || (this.mouseButton === 3 && ev.ctrlKey)) { // Slab
                this.slabNear = this.cslabNear + dx * 100;
                this.slabFar = this.cslabFar + dy * 100;
            } else if (mode === 2 || this.mouseButton === 3 || ev.shiftKey) { // Zoom
                scaleFactor = (this.rotationGroup.position.z - this.CAMERA_Z) * 0.85;
                if (scaleFactor < 80) { scaleFactor = 80; }
                this.rotationGroup.position.z = this.cz - dy * scaleFactor;
            } else if (mode === 1 || this.mouseButton === 2 || ev.ctrlKey) { // Translate
                scaleFactor = (this.rotationGroup.position.z - this.CAMERA_Z) * 0.85;
                if (scaleFactor < 20) { scaleFactor = 20; }
                translationByScreen = new TV3(-dx * scaleFactor, -dy * scaleFactor, 0);
                q = this.rotationGroup.quaternion;
                qinv = new THREE.Quaternion(q.x, q.y, q.z, q.w).inverse().normalize();
                translation = qinv.multiplyVector3(translationByScreen);
                this.modelGroup.position.x = this.currentModelPos.x + translation.x;
                this.modelGroup.position.y = this.currentModelPos.y + translation.y;
                this.modelGroup.position.z = this.currentModelPos.z + translation.z;
            } else if ((mode === 0 || this.mouseButton === 1) && r !== 0) { // Rotate
                rs = Math.sin(r * Math.PI) / r;
                this.dq.x = Math.cos(r * Math.PI);
                this.dq.y = 0;
                this.dq.z =  rs * dx;
                this.dq.w =  rs * dy;
                this.rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0);
                this.rotationGroup.quaternion.multiplySelf(this.dq);
                this.rotationGroup.quaternion.multiplySelf(this.cq);
            }
            this.show();
        }.bind(this));
    };


    GLmol.prototype.show = function () {
        if (!this.scene) {
            return;
        }

        console.time("rendered");
        this.setSlabAndFog();
        this.renderer.render(this.scene, this.camera);
        console.timeEnd("rendered");
    };

    // For scripting
    GLmol.prototype.doFunc = function (func) {
        func(this);
    };

    // Expose GLmol to the global object
    window.GLmol = GLmol;


}(window));
