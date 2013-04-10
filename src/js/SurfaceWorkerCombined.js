//dkoes - calls protein surface as a webworker, returns the faces and vertices

/*
 * These tables are based off those by Paul Bourke and Geoffrey Heller:
 * http://paulbourke.net/geometry/polygonise/
 * http://paulbourke.net/geometry/polygonise/table2.txt
 * 
 * However, they have been substantially modified to reflect a more 
 * sensible corner numbering scheme and the discrete nature of our voxel data
 * (resulting in fewer faces).
 */
var MarchingCube = MarchingCube || {};

MarchingCube.edgeTable = new Uint32Array([ 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0xb00, 0x0, 0x0, 0x0, 0x700, 0x0, 0xd00, 0xe00, 0xf00, 0x0, 0x0, 0x0,
		0x8a, 0x0, 0x15, 0x0, 0x86, 0x0, 0x0, 0x0, 0x28c, 0x0, 0x813, 0xf19,
		0xe10, 0x0, 0x0, 0x0, 0x2a, 0x0, 0x0, 0x0, 0x126, 0x0, 0x0, 0x15, 0x1c,
		0x0, 0xf23, 0x419, 0xd20, 0x0, 0xa8, 0xa2, 0xaa, 0x0, 0x285, 0x9ab,
		0x8a2, 0x0, 0x2af, 0x125, 0xac, 0xfaa, 0xea3, 0xda9, 0xca0, 0x0, 0x0,
		0x0, 0x0, 0x0, 0x45, 0x0, 0x384, 0x0, 0x0, 0x0, 0x700, 0x8a, 0x83,
		0x648, 0x780, 0x0, 0x51, 0x0, 0x81a, 0x54, 0x55, 0x54, 0x56, 0x0, 0x51,
		0x0, 0xe5c, 0x14a, 0x451, 0x759, 0x650, 0x0, 0x0, 0x0, 0x2a, 0x0, 0x45,
		0x0, 0x1f6, 0x0, 0x0, 0x15, 0xdfc, 0x8a, 0x7f3, 0x4f9, 0x5f0, 0xb00,
		0x68, 0x921, 0x6a, 0x348, 0x245, 0x16f, 0x66, 0xb00, 0xe6f, 0xd65,
		0xc6c, 0x76a, 0x663, 0x569, 0x460, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
		0xf46, 0x0, 0x0, 0x45, 0x24c, 0x2a, 0x823, 0x29, 0xb40, 0x0, 0x0, 0x0,
		0x6ba, 0x0, 0x8f5, 0xfff, 0xef6, 0x0, 0xff, 0x2f5, 0x2fc, 0x9ea, 0x8f3,
		0xbf9, 0xaf0, 0x0, 0x0, 0x51, 0x152, 0x0, 0xf55, 0x45f, 0xd56, 0x54,
		0x357, 0x55, 0x154, 0x852, 0xb53, 0x59, 0x950, 0x700, 0x2c8, 0xc2,
		0x48a, 0xfc4, 0xec5, 0xdcf, 0xcc6, 0x2c4, 0x2cf, 0xc5, 0xcc, 0xbca,
		0xac3, 0x9c9, 0x8c0, 0x0, 0x0, 0x0, 0x0, 0xa8, 0x1a4, 0xa8, 0x7a6,
		0xa2, 0xa2, 0x2a4, 0xbac, 0xaa, 0xa3, 0x2a8, 0x3a0, 0xd00, 0xc18,
		0xd00, 0xe3a, 0x34, 0x35, 0x73f, 0x636, 0x924, 0x83f, 0xb35, 0xa3c,
		0x12a, 0x33, 0x339, 0x230, 0xe00, 0xe00, 0xc12, 0xd9a, 0x684, 0x795,
		0x49f, 0x596, 0x92, 0xb9f, 0x815, 0x99c, 0x9a, 0x393, 0x99, 0x190,
		0xf00, 0xe08, 0xd01, 0xc0a, 0x704, 0x605, 0x50f, 0x406, 0xb02, 0xa0f,
		0x905, 0x80c, 0x30a, 0x203, 0x109, 0x0 ]);

MarchingCube.triTable = [ [], [], [], [], [], [], [], [ 11, 9, 8 ], [], [], [],
		[ 8, 10, 9 ], [], [ 10, 8, 11 ], [ 9, 11, 10 ],
		[ 8, 10, 9, 8, 11, 10 ], [], [], [], [ 1, 7, 3 ], [], [ 4, 2, 0 ], [],
		[ 2, 1, 7 ], [], [], [], [ 2, 7, 3, 2, 9, 7 ], [],
		[ 1, 4, 11, 1, 0, 4 ], [ 3, 8, 0, 11, 9, 4, 11, 10, 9 ],
		[ 4, 11, 9, 11, 10, 9 ], [], [], [], [ 5, 3, 1 ], [], [], [],
		[ 2, 5, 8, 2, 1, 5 ], [], [], [ 2, 4, 0 ], [ 3, 2, 4 ], [],
		[ 0, 9, 1, 8, 10, 5, 8, 11, 10 ], [ 3, 4, 0, 3, 10, 4 ],
		[ 5, 8, 10, 8, 11, 10 ], [], [ 3, 5, 7 ], [ 7, 1, 5 ],
		[ 1, 7, 3, 1, 5, 7 ], [], [ 9, 2, 0, 9, 7, 2 ],
		[ 0, 3, 8, 1, 7, 11, 1, 5, 7 ], [ 11, 1, 7, 1, 5, 7 ], [],
		[ 9, 1, 0, 5, 3, 2, 5, 7, 3 ], [ 8, 2, 5, 8, 0, 2 ],
		[ 2, 5, 3, 5, 7, 3 ], [ 3, 9, 1, 3, 8, 9, 7, 11, 10, 7, 10, 5 ],
		[ 9, 1, 0, 10, 7, 11, 10, 5, 7 ], [ 3, 8, 0, 7, 10, 5, 7, 11, 10 ],
		[ 11, 5, 7, 11, 10, 5 ], [], [], [], [], [], [ 0, 6, 2 ], [],
		[ 7, 2, 9, 7, 9, 8 ], [], [], [], [ 8, 10, 9 ], [ 7, 1, 3 ],
		[ 7, 1, 0 ], [ 6, 9, 3, 6, 10, 9 ], [ 7, 10, 8, 10, 9, 8 ], [],
		[ 6, 0, 4 ], [], [ 11, 1, 4, 11, 3, 1 ], [ 2, 4, 6 ],
		[ 2, 0, 4, 2, 4, 6 ], [ 2, 4, 6 ], [ 1, 4, 2, 4, 6, 2 ], [],
		[ 6, 0, 4 ], [], [ 2, 11, 3, 6, 9, 4, 6, 10, 9 ], [ 8, 6, 1, 8, 1, 3 ],
		[ 10, 0, 6, 0, 4, 6 ], [ 8, 0, 3, 9, 6, 10, 9, 4, 6 ],
		[ 10, 4, 6, 10, 9, 4 ], [], [], [], [ 5, 3, 1 ], [], [ 0, 6, 2 ], [],
		[ 7, 4, 8, 5, 2, 1, 5, 6, 2 ], [], [], [ 2, 4, 0 ],
		[ 7, 4, 8, 2, 11, 3, 10, 5, 6 ], [ 7, 1, 3 ],
		[ 5, 6, 10, 0, 9, 1, 8, 7, 4 ], [ 5, 6, 10, 7, 0, 3, 7, 4, 0 ],
		[ 10, 5, 6, 4, 8, 7 ], [ 9, 11, 8 ], [ 3, 5, 6 ],
		[ 0, 5, 11, 0, 11, 8 ], [ 6, 3, 5, 3, 1, 5 ], [ 3, 9, 6, 3, 8, 9 ],
		[ 9, 6, 0, 6, 2, 0 ], [ 0, 3, 8, 2, 5, 6, 2, 1, 5 ],
		[ 1, 6, 2, 1, 5, 6 ], [ 9, 11, 8 ], [ 1, 0, 9, 6, 10, 5, 11, 3, 2 ],
		[ 6, 10, 5, 2, 8, 0, 2, 11, 8 ], [ 3, 2, 11, 10, 5, 6 ],
		[ 10, 5, 6, 9, 3, 8, 9, 1, 3 ], [ 0, 9, 1, 5, 6, 10 ],
		[ 8, 0, 3, 10, 5, 6 ], [ 10, 5, 6 ], [], [], [], [], [], [], [],
		[ 1, 10, 2, 9, 11, 6, 9, 8, 11 ], [], [], [ 6, 0, 2 ],
		[ 3, 6, 9, 3, 2, 6 ], [ 3, 5, 1 ], [ 0, 5, 1, 0, 11, 5 ], [ 0, 3, 5 ],
		[ 6, 9, 11, 9, 8, 11 ], [], [], [], [ 4, 5, 9, 7, 1, 10, 7, 3, 1 ], [],
		[ 11, 6, 7, 2, 4, 5, 2, 0, 4 ],
		[ 11, 6, 7, 8, 0, 3, 1, 10, 2, 9, 4, 5 ],
		[ 6, 7, 11, 1, 10, 2, 9, 4, 5 ], [],
		[ 4, 1, 0, 4, 5, 1, 6, 7, 3, 6, 3, 2 ], [ 9, 4, 5, 0, 6, 7, 0, 2, 6 ],
		[ 4, 5, 9, 6, 3, 2, 6, 7, 3 ], [ 6, 7, 11, 5, 3, 8, 5, 1, 3 ],
		[ 6, 7, 11, 4, 1, 0, 4, 5, 1 ], [ 4, 5, 9, 3, 8, 0, 11, 6, 7 ],
		[ 9, 4, 5, 7, 11, 6 ], [], [], [ 0, 6, 4 ], [ 8, 6, 4, 8, 1, 6 ], [],
		[ 0, 10, 2, 0, 9, 10, 4, 8, 11, 4, 11, 6 ],
		[ 10, 2, 1, 6, 0, 3, 6, 4, 0 ], [ 10, 2, 1, 11, 4, 8, 11, 6, 4 ],
		[ 4, 2, 6 ], [ 1, 0, 9, 2, 4, 8, 2, 6, 4 ], [ 2, 4, 0, 2, 6, 4 ],
		[ 8, 2, 4, 2, 6, 4 ], [ 11, 4, 1, 11, 6, 4 ],
		[ 0, 9, 1, 4, 11, 6, 4, 8, 11 ], [ 3, 6, 0, 6, 4, 0 ],
		[ 8, 6, 4, 8, 11, 6 ], [ 10, 8, 9 ], [ 6, 3, 9, 6, 7, 3 ], [ 6, 7, 1 ],
		[ 10, 7, 1, 7, 3, 1 ], [ 7, 11, 6, 8, 10, 2, 8, 9, 10 ],
		[ 11, 6, 7, 10, 0, 9, 10, 2, 0 ], [ 2, 1, 10, 7, 11, 6, 8, 0, 3 ],
		[ 1, 10, 2, 6, 7, 11 ], [ 7, 2, 6, 7, 9, 2 ],
		[ 1, 0, 9, 3, 6, 7, 3, 2, 6 ], [ 7, 0, 6, 0, 2, 6 ],
		[ 2, 7, 3, 2, 6, 7 ], [ 7, 11, 6, 3, 9, 1, 3, 8, 9 ],
		[ 9, 1, 0, 11, 6, 7 ], [ 0, 3, 8, 11, 6, 7 ], [ 11, 6, 7 ], [], [], [],
		[], [ 5, 3, 7 ], [ 8, 5, 2, 8, 7, 5 ], [ 5, 3, 7 ],
		[ 1, 10, 2, 5, 8, 7, 5, 9, 8 ], [ 1, 7, 5 ], [ 1, 7, 5 ],
		[ 9, 2, 7, 9, 7, 5 ], [ 11, 3, 2, 8, 5, 9, 8, 7, 5 ],
		[ 1, 3, 7, 1, 7, 5 ], [ 0, 7, 1, 7, 5, 1 ], [ 9, 3, 5, 3, 7, 5 ],
		[ 9, 7, 5, 9, 8, 7 ], [ 8, 10, 11 ], [ 3, 4, 10, 3, 10, 11 ],
		[ 8, 10, 11 ], [ 5, 9, 4, 1, 11, 3, 1, 10, 11 ], [ 2, 4, 5 ],
		[ 5, 2, 4, 2, 0, 4 ], [ 0, 3, 8, 5, 9, 4, 10, 2, 1 ],
		[ 2, 1, 10, 9, 4, 5 ], [ 2, 8, 5, 2, 11, 8 ],
		[ 3, 2, 11, 1, 4, 5, 1, 0, 4 ], [ 9, 4, 5, 8, 2, 11, 8, 0, 2 ],
		[ 11, 3, 2, 9, 4, 5 ], [ 8, 5, 3, 5, 1, 3 ], [ 5, 0, 4, 5, 1, 0 ],
		[ 3, 8, 0, 4, 5, 9 ], [ 9, 4, 5 ], [ 11, 9, 10 ], [ 11, 9, 10 ],
		[ 1, 11, 4, 1, 10, 11 ], [ 8, 7, 4, 11, 1, 10, 11, 3, 1 ],
		[ 2, 7, 9, 2, 9, 10 ], [ 4, 8, 7, 0, 10, 2, 0, 9, 10 ],
		[ 2, 1, 10, 0, 7, 4, 0, 3, 7 ], [ 10, 2, 1, 8, 7, 4 ], [ 1, 7, 4 ],
		[ 3, 2, 11, 4, 8, 7, 9, 1, 0 ], [ 11, 4, 2, 4, 0, 2 ],
		[ 2, 11, 3, 7, 4, 8 ], [ 4, 1, 7, 1, 3, 7 ], [ 1, 0, 9, 8, 7, 4 ],
		[ 3, 4, 0, 3, 7, 4 ], [ 8, 7, 4 ], [ 8, 9, 10, 8, 10, 11 ],
		[ 3, 9, 11, 9, 10, 11 ], [ 0, 10, 8, 10, 11, 8 ],
		[ 10, 3, 1, 10, 11, 3 ], [ 2, 8, 10, 8, 9, 10 ], [ 9, 2, 0, 9, 10, 2 ],
		[ 8, 0, 3, 1, 10, 2 ], [ 10, 2, 1 ], [ 1, 11, 9, 11, 8, 9 ],
		[ 11, 3, 2, 0, 9, 1 ], [ 11, 0, 2, 11, 8, 0 ], [ 11, 3, 2 ],
		[ 8, 1, 3, 8, 9, 1 ], [ 9, 1, 0 ], [ 8, 0, 3 ], [] ];

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

=======

TODO: Improved performance on Firefox
      Reduce memory consumption
      Refactor!
 */

// dkoes
// Surface calculations.  This must be safe to use within a web worker.
if (typeof console === 'undefined') {
	// this should only be true inside of a webworker
	console = {
		log : function() {
		}
	};
}

var ProteinSurface = (function() {
	// constants for vpbits bitmasks
	var INOUT = 1;
	var ISDONE = 2;
	var ISBOUND = 4;

	var ptranx = 0, ptrany = 0, ptranz = 0;
	var probeRadius = 1.4;
	var scaleFactor = 2; // 2 is .5A grid; if this is made user configurable,
							// also have to adjust offset used to find non-shown
							// atoms
	var pHeight = 0, pWidth = 0, pLength = 0;
	var cutRadius = 0;
	var vpBits = null; // uint8 array of bitmasks
	var vpDistance = null; // floatarray
	var vpAtomID = null; // intarray
	var vertnumber = 0, facenumber = 0;
	var pminx = 0, pminy = 0, pminz = 0, pmaxx = 0, pmaxy = 0, pmaxz = 0;

	var vdwRadii = {
			"H" : 1.2,
			"Li" : 1.82,
			"Na" : 2.27,
			"K" : 2.75,
			"C" : 1.7,
			"N" : 1.55,
			"O" : 1.52,
			"F" : 1.47,
			"P" : 1.80,
			"S" : 1.80,
			"CL" : 1.75,
			"BR" : 1.85,
			"SE" : 1.90,
			"ZN" : 1.39,
			"CU" : 1.4,
			"NI" : 1.63,
			"X" : 2
		};

	var getVDWIndex = function(atom) {
		if(!atom.elem || typeof(vdwRadii[atom.elem]) == "undefined") {
			return "X";
		}
		return atom.elem;
	};
	
	var depty = {}, widxz = {};
	var faces, verts;
	var nb = [ [ 1, 0, 0 ], [ -1, 0, 0 ], [ 0, 1, 0 ], [ 0, -1, 0 ],
			[ 0, 0, 1 ], [ 0, 0, -1 ], [ 1, 1, 0 ], [ 1, -1, 0 ], [ -1, 1, 0 ],
			[ -1, -1, 0 ], [ 1, 0, 1 ], [ 1, 0, -1 ], [ -1, 0, 1 ],
			[ -1, 0, -1 ], [ 0, 1, 1 ], [ 0, 1, -1 ], [ 0, -1, 1 ],
			[ 0, -1, -1 ], [ 1, 1, 1 ], [ 1, 1, -1 ], [ 1, -1, 1 ],
			[ -1, 1, 1 ], [ 1, -1, -1 ], [ -1, -1, 1 ], [ -1, 1, -1 ],
			[ -1, -1, -1 ] ];

	var origextent;

	var inOrigExtent = function(x, y, z) {
		if (x < origextent[0][0] || x > origextent[1][0])
			return false;
		if (y < origextent[0][1] || y > origextent[1][1])
			return false;
		if (z < origextent[0][2] || z > origextent[1][2])
			return false;
		return true;
	};

	this.getFacesAndVertices = function(atomlist) {
		var atomsToShow = new Object();
		for ( var i = 0, lim = atomlist.length; i < lim; i++)
			atomsToShow[atomlist[i]] = true;
		var vertices = verts;
		for (i = 0; i < vertices.length; i++) {
			vertices[i].x = vertices[i].x / scaleFactor - ptranx;
			vertices[i].y = vertices[i].y / scaleFactor - ptrany;
			vertices[i].z = vertices[i].z / scaleFactor - ptranz;
		}

		var finalfaces = []
		for ( var i = 0; i < faces.length; i++) {
			var f = faces[i];
			var a = vertices[f.a].atomid, b = vertices[f.b].atomid, c = vertices[f.c].atomid;

			// must be a unique face for each atom
			var which = a;

			if (b < which)
				which = b;
			if (c < which)
				which = c;
			if (!atomsToShow[which]) {
				continue;
			}
			var av = vertices[f.a];
			var bv = vertices[f.b];
			var cv = vertices[f.c];

			if (f.a != f.bb && f.b != f.c && f.a != f.c)
				finalfaces.push(f);
		}

		//try to help the garbage collector
		vpBits = null; // uint8 array of bitmasks
		vpDistance = null; // floatarray
		vpAtomID = null; // intarray
		
		return {
			vertices : vertices,
			faces : finalfaces
		};
	};

	this.laplaciansmooth = function(numiter) {
		var tps = new Array(verts.length);
		for ( var i = 0; i < verts.length; i++)
			tps[i] = {
				x : 0,
				y : 0,
				z : 0
			};
		var vertdeg = new Array(20);
		var flagvert;
		for ( var i = 0; i < 20; i++)
			vertdeg[i] = new Array(verts.length);
		for ( var i = 0; i < verts.length; i++)
			vertdeg[0][i] = 0;
		for ( var i = 0; i < faces.length; i++) {
			flagvert = true;
			for ( var j = 0; j < vertdeg[0][faces[i].a]; j++) {
				if (faces[i].b == vertdeg[j + 1][faces[i].a]) {
					flagvert = false;
					break;
				}
			}
			if (flagvert) {
				vertdeg[0][faces[i].a]++;
				vertdeg[vertdeg[0][faces[i].a]][faces[i].a] = faces[i].b;
			}
			flagvert = true;
			for ( var j = 0; j < vertdeg[0][faces[i].a]; j++) {
				if (faces[i].c == vertdeg[j + 1][faces[i].a]) {
					flagvert = false;
					break;
				}
			}
			if (flagvert) {
				vertdeg[0][faces[i].a]++;
				vertdeg[vertdeg[0][faces[i].a]][faces[i].a] = faces[i].c;
			}
			// b
			flagvert = true;
			for (j = 0; j < vertdeg[0][faces[i].b]; j++) {
				if (faces[i].a == vertdeg[j + 1][faces[i].b]) {
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
				if (faces[i].c == vertdeg[j + 1][faces[i].b]) {
					flagvert = false;
					break;
				}
			}
			if (flagvert) {
				vertdeg[0][faces[i].b]++;
				vertdeg[vertdeg[0][faces[i].b]][faces[i].b] = faces[i].c;
			}
			// c
			flagvert = true;
			for (j = 0; j < vertdeg[0][faces[i].c]; j++) {
				if (faces[i].a == vertdeg[j + 1][faces[i].c]) {
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
				if (faces[i].b == vertdeg[j + 1][faces[i].c]) {
					flagvert = false;
					break;
				}
			}
			if (flagvert) {
				vertdeg[0][faces[i].c]++;
				vertdeg[vertdeg[0][faces[i].c]][faces[i].c] = faces[i].b;
			}
		}

		var wt = 1.00;
		var wt2 = 0.50;
		var ssign;
		var outwt = 0.75 / (scaleFactor + 3.5); // area-preserving
		for ( var k = 0; k < numiter; k++) {
			for ( var i = 0; i < verts.length; i++) {
				if (vertdeg[0][i] < 3) {
					tps[i].x = verts[i].x;
					tps[i].y = verts[i].y;
					tps[i].z = verts[i].z;
				} else if (vertdeg[0][i] == 3 || vertdeg[0][i] == 4) {
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
					for ( var j = 0; j < vertdeg[0][i]; j++) {
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
			for ( var i = 0; i < verts.length; i++) {
				verts[i].x = tps[i].x;
				verts[i].y = tps[i].y;
				verts[i].z = tps[i].z;
			}
			/*
			 * computenorm(); for (var i = 0; i < vertnumber; i++) { if
			 * (verts[i].inout) ssign = 1; else ssign = -1; verts[i].x += ssign *
			 * outwt * verts[i].pn.x; verts[i].y += ssign * outwt *
			 * verts[i].pn.y; verts[i].z += ssign * outwt * verts[i].pn.z; }
			 */
		}
	};

	this.initparm = function(extent, btype) {
		var margin = (1 / scaleFactor) * 5.5; // need margine to avoid
												// boundary/round off effects
		origextent = extent;
		pminx = extent[0][0], pmaxx = extent[1][0];
		pminy = extent[0][1], pmaxy = extent[1][1];
		pminz = extent[0][2], pmaxz = extent[1][2];

		if (!btype) {
			pminx -= margin;
			pminy -= margin;
			pminz -= margin;
			pmaxx += margin;
			pmaxy += margin;
			pmaxz += margin;
		} else {
			pminx -= probeRadius + margin;
			pminy -= probeRadius + margin;
			pminz -= probeRadius + margin;
			pmaxx += probeRadius + margin;
			pmaxy += probeRadius + margin;
			pmaxz += probeRadius + margin;
		}

		pminx = Math.floor(pminx * scaleFactor) / scaleFactor;
		pminy = Math.floor(pminy * scaleFactor) / scaleFactor;
		pminz = Math.floor(pminz * scaleFactor) / scaleFactor;
		pmaxx = Math.ceil(pmaxx * scaleFactor) / scaleFactor;
		pmaxy = Math.ceil(pmaxy * scaleFactor) / scaleFactor;
		pmaxz = Math.ceil(pmaxz * scaleFactor) / scaleFactor;

		ptranx = -pminx;
		ptrany = -pminy;
		ptranz = -pminz;

		pLength = Math.ceil(scaleFactor * (pmaxx - pminx)) + 1;
		pWidth = Math.ceil(scaleFactor * (pmaxy - pminy)) + 1;
		pHeight = Math.ceil(scaleFactor * (pmaxz - pminz)) + 1;

		this.boundingatom(btype);
		cutRadis = probeRadius * scaleFactor;

		vpBits = new Uint8Array(pLength * pWidth * pHeight);
		vpDistance = new Float64Array(pLength * pWidth * pHeight); // float 32
		// doesn't
		// play
		// nicely
		// with
		// native
		// floats
		vpAtomID = new Int32Array(pLength * pWidth * pHeight);
		console.log("Box size: ", pLength, pWidth, pHeight, vpBits.length);
	};

	this.boundingatom = function(btype) {
		var tradius = [];
		var txz, tdept, sradius, idx;
		flagradius = btype;

		for ( var i in vdwRadii) {
			if(!vdwRadii.hasOwnProperty(i))
				continue;
			var r = vdwRadii[i];
			if (!btype)
				tradius[i] = r * scaleFactor + 0.5;
			else
				tradius[i] = (r + probeRadius) * scaleFactor + 0.5;

			sradius = tradius[i] * tradius[i];
			widxz[i] = Math.floor(tradius[i]) + 1;
			depty[i] = new Array(widxz[i] * widxz[i]);
			indx = 0;
			for (j = 0; j < widxz[i]; j++) {
				for (k = 0; k < widxz[i]; k++) {
					txz = j * j + k * k;
					if (txz > sradius)
						depty[i][indx] = -1; // outside
					else {
						tdept = Math.sqrt(sradius - txz);
						depty[i][indx] = Math.floor(tdept);
					}
					indx++;
				}
			}
		}
	};

	this.fillvoxels = function(atoms, atomlist) { // (int seqinit,int
		// seqterm,bool
		// atomtype,atom*
		// proseq,bool bcolor)
		var i;
		for ( var i = 0, lim = vpBits.length; i < lim; i++) {
			vpBits[i] = 0;
			vpDistance[i] = -1.0;
			vpAtomID[i] = -1;
		}

		for (i in atomlist) {
			atom = atoms[atomlist[i]];
			if (atom == undefined)
				continue;
			this.fillAtom(atom, atoms);
		}

		for (i = 0, lim = vpBits.length; i < lim; i++)
			if (vpBits[i] & INOUT)
				vpBits[i] |= ISDONE;

	};


	this.fillAtom = function(atom, atoms) {
		var cx, cy, cz, ox, oy, oz;
		cx = Math.floor(0.5 + scaleFactor * (atom.x + ptranx));
		cy = Math.floor(0.5 + scaleFactor * (atom.y + ptrany));
		cz = Math.floor(0.5 + scaleFactor * (atom.z + ptranz));

		var at = getVDWIndex(atom);
		var nind = 0;
		var cnt = 0;

		for (i = 0, n = widxz[at]; i < n; i++) {
			for (j = 0; j < n; j++) {
				if (depty[at][nind] != -1) {
					for (ii = -1; ii < 2; ii++) {
						for (jj = -1; jj < 2; jj++) {
							for (kk = -1; kk < 2; kk++) {
								if (ii != 0 && jj != 0 && kk != 0) {
									mi = ii * i;
									mk = kk * j;
									for (k = 0; k <= depty[at][nind]; k++) {
										mj = k * jj;
										si = cx + mi;
										sj = cy + mj;
										sk = cz + mk;
										if (si < 0 || sj < 0 || sk < 0
												|| si >= pLength
												|| sj >= pWidth
												|| sk >= pHeight)
											continue;
										var index = si * pWidth * pHeight + sj
												* pHeight + sk;

										if (!(vpBits[index] & INOUT)) {
											vpBits[index] |= INOUT;
											vpAtomID[index] = atom.serial;
										} else {
											var atom2 = atoms[vpAtomID[index]];
											ox = Math.floor(0.5 + scaleFactor
													* (atom2.x + ptranx));
											oy = Math.floor(0.5 + scaleFactor
													* (atom2.y + ptrany));
											oz = Math.floor(0.5 + scaleFactor
													* (atom2.z + ptranz));
											if (mi * mi + mj * mj + mk * mk < ox
													* ox + oy * oy + oz * oz)
												vpAtomID[index] = atom.serial;
										}

									}// k
								}// if
							}// kk
						}// jj
					}// ii
				}// if
				nind++;
			}// j
		}// i
	};

	this.fillvoxelswaals = function(atoms, atomlist) {
		for ( var i = 0, lim = vpBits.length; i < lim; i++)
			vpBits[i] &= ~ISDONE; // not isdone

		for (i in atomlist) {
			atom = atoms[atomlist[i]];
			if (atom == undefined)
				continue;

			this.fillAtomWaals(atom, atoms);
		}
	};

	this.fillAtomWaals = function(atom, atoms) {
		var cx, cy, cz, ox, oy, oz, nind = 0;
		cx = Math.floor(0.5 + scaleFactor * (atom.x + ptranx));
		cy = Math.floor(0.5 + scaleFactor * (atom.y + ptrany));
		cz = Math.floor(0.5 + scaleFactor * (atom.z + ptranz));

		var at = getVDWIndex(atom);
		var i, j,n;
		for (i = 0, n = widxz[at]; i < n; i++) {
			for (j = 0; j < n; j++) {
				if (depty[at][nind] != -1) {
					for (ii = -1; ii < 2; ii++) {
						for (jj = -1; jj < 2; jj++) {
							for (kk = -1; kk < 2; kk++) {
								if (ii != 0 && jj != 0 && kk != 0) {
									mi = ii * i;
									mk = kk * j;
									for (k = 0; k <= depty[at][nind]; k++) {
										mj = k * jj;
										si = cx + mi;
										sj = cy + mj;
										sk = cz + mk;
										if (si < 0 || sj < 0 || sk < 0
												|| si >= pLength
												|| sj >= pWidth
												|| sk >= pHeight)
											continue;
										var index = si * pWidth * pHeight + sj
												* pHeight + sk;
										if (!(vpBits[index] & ISDONE)) {
											vpBits[index] |= ISDONE;
											vpAtomID[index] = atom.serial;
										} else {
											var atom2 = atoms[vpAtomID[index]];
											ox = Math.floor(0.5 + scaleFactor
													* (atom2.x + ptranx));
											oy = Math.floor(0.5 + scaleFactor
													* (atom2.y + ptrany));
											oz = Math.floor(0.5 + scaleFactor
													* (atom2.z + ptranz));
											if (mi * mi + mj * mj + mk * mk < ox
													* ox + oy * oy + oz * oz)
												vpAtomID[index] = atom.serial;
										}

									}// k
								}// if
							}// kk
						}// jj
					}// ii
				}// if
				nind++;
			}// j
		}// i
	};

	this.buildboundary = function() {
		for (i = 0; i < pLength; i++) {
			for (j = 0; j < pHeight; j++) {
				for (k = 0; k < pWidth; k++) {
					var index = i * pWidth * pHeight + k * pHeight + j;
					if (vpBits[index] & INOUT) {
						var flagbound = false;
						var ii = 0;
						while (!flagbound && ii < 26) {
							var ti = i + nb[ii][0], tj = j + nb[ii][2], tk = k
									+ nb[ii][1];
							if (ti > -1
									&& ti < pLength
									&& tk > -1
									&& tk < pWidth
									&& tj > -1
									&& tj < pHeight
									&& !(vpBits[ti * pWidth * pHeight + tk
											* pHeight + tj] & INOUT)) {
								vpBits[index] |= ISBOUND;
								flagbound = true;
							} else
								ii++;
						}
					}
				}
			}
		}
	};

	// a little class for 3d array, should really generalize this and
	// use throughout...
	var PointGrid = function(length, width, height) {
		// the standard says this is zero initialized
		var data = new Int32Array(length * width * height * 3);

		// set position x,y,z to pt, which has ix,iy,and iz
		this.set = function(x, y, z, pt) {
			var index = ((((x * width) + y) * height) + z) * 3;
			data[index] = pt.ix;
			data[index + 1] = pt.iy;
			data[index + 2] = pt.iz;
		};

		// return point at x,y,z
		this.get = function(x, y, z) {
			var index = ((((x * width) + y) * height) + z) * 3;
			return {
				ix : data[index],
				iy : data[index + 1],
				iz : data[index + 2]
			};
		};
	};

	this.fastdistancemap = function() {
		var positin = 0, positout = 0, eliminate = 0;
		var certificate;
		var i, j, k;
		totalsurfacevox = 0;
		totalinnervox = 0;

		var boundPoint = new PointGrid(pLength, pWidth, pHeight);

		for (i = 0; i < pLength; i++) {
			for (j = 0; j < pWidth; j++) {
				for (k = 0; k < pHeight; k++) {
					var index = i * pWidth * pHeight + j * pHeight + k;
					vpBits[index] &= ~ISDONE; // isdone = false
					if (vpBits[index] & INOUT) {
						if (vpBits[index] & ISBOUND) {
							totalsurfacevox++;
							boundPoint.set(i, j, k, {
								ix : i,
								iy : j,
								iz : k
							});
							vpDistance[index] = 0;
							vpBits[index] |= ISDONE;
						} else {
							totalinnervox++;
						}
					}
				}
			}
		}

		inarray = new Array();
		outarray = new Array();
		var positin = 0, positout = 0;

		for (i = 0; i < pLength; i++) {
			for (j = 0; j < pWidth; j++) {
				for (k = 0; k < pHeight; k++) {
					var index = i * pWidth * pHeight + j * pHeight + k;
					if (vpBits[index] & ISBOUND) {
						inarray.push({
							ix : i,
							iy : j,
							iz : k
						});
						positin++;
						vpBits[index] &= ~ISBOUND;
					}
				}
			}
		}

		do {
			positout = this.fastoneshell(positin, boundPoint);
			positin = 0;
			inarray = [];
			for (i = 0; i < positout; i++) {
				var index = pWidth * pHeight * outarray[i].ix + pHeight
						* outarray[i].iy + outarray[i].iz;
				vpBits[index] &= ~ISBOUND;
				if (vpDistance[index] <= 1.02 * cutRadis) {
					inarray.push({
						ix : outarray[i].ix,
						iy : outarray[i].iy,
						iz : outarray[i].iz
					});
					// inarray[positin].ix=outarray[i].ix;
					// inarray[positin].iy=outarray[i].iy;
					// inarray[positin].iz=outarray[i].iz;
					positin++;
				}
			}
		} while (positin != 0);

		var cutsf = scaleFactor - 0.5;
		if (cutsf < 0)
			cutsf = 0;
		for (i = 0; i < pLength; i++) {
			for (j = 0; j < pWidth; j++) {
				for (k = 0; k < pHeight; k++) {
					var index = i * pWidth * pHeight + j * pHeight + k;
					vpBits[index] &= ~ISBOUND;
					// ses solid
					if (vpBits[index] & INOUT) {
						if (!(vpBits[index] & ISDONE)
								|| ((vpBits[index] & ISDONE) && vpDistance[index] >= cutRadis
										- 0.50 / (0.1 + cutsf))) {
							vpBits[index] |= ISBOUND;
						}
					}
				}
			}
		}
		inarray = [];
		outarray = [];
		boundPoint = null;
	};

	this.fastoneshell = function(number, boundPoint) { // (int* innum,int
		// *allocout,voxel2
		// ***boundPoint, int*
		// outnum, int *elimi)
		var positout = 0;
		var tx, ty, tz;
		var dx, dy, dz;
		var square;
		if (number == 0)
			return 0;
		outarray = [];

		tnv = {
			ix : -1,
			iy : -1,
			iz : -1
		};
		for ( var i = 0; i < number; i++) {
			tx = inarray[i].ix;
			ty = inarray[i].iy;
			tz = inarray[i].iz;
			var bp = boundPoint.get(tx, ty, tz);

			for ( var j = 0; j < 6; j++) {
				tnv.ix = tx + nb[j][0];
				tnv.iy = ty + nb[j][1];
				tnv.iz = tz + nb[j][2];
				var index = tnv.ix * pWidth * pHeight + pHeight * tnv.iy
						+ tnv.iz;
				if (tnv.ix < pLength && tnv.ix > -1 && tnv.iy < pWidth
						&& tnv.iy > -1 && tnv.iz < pHeight && tnv.iz > -1
						&& (vpBits[index] & INOUT) && !(vpBits[index] & ISDONE)) {

					boundPoint.set(tnv.ix, tnv.iy, tz + nb[j][2], bp);
					dx = tnv.ix - bp.ix;
					dy = tnv.iy - bp.iy;
					dz = tnv.iz - bp.iz;
					var square = dx * dx + dy * dy + dz * dz;
					vpDistance[index] = Math.sqrt(square);
					vpBits[index] |= ISDONE;
					vpBits[index] |= ISBOUND;

					outarray.push({
						ix : tnv.ix,
						iy : tnv.iy,
						iz : tnv.iz
					});
					positout++;
				} else if (tnv.ix < pLength && tnv.ix > -1 && tnv.iy < pWidth
						&& tnv.iy > -1 && tnv.iz < pHeight && tnv.iz > -1
						&& (vpBits[index] & INOUT) && (vpBits[index] & ISDONE)) {

					dx = tnv.ix - bp.ix;
					dy = tnv.iy - bp.iy;
					dz = tnv.iz - bp.iz;
					square = dx * dx + dy * dy + dz * dz;
					square = Math.sqrt(square);
					if (square < vpDistance[index]) {
						boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);

						vpDistance[index] = square;
						if (!(vpBits[index] & ISBOUND)) {
							vpBits[index] |= ISBOUND;
							outarray.push({
								ix : tnv.ix,
								iy : tnv.iy,
								iz : tnv.iz
							});
							positout++;
						}
					}
				}
			}
		}

		// console.log("part1", positout);

		for (i = 0; i < number; i++) {
			tx = inarray[i].ix;
			ty = inarray[i].iy;
			tz = inarray[i].iz;
			var bp = boundPoint.get(tx, ty, tz);

			for (j = 6; j < 18; j++) {
				tnv.ix = tx + nb[j][0];
				tnv.iy = ty + nb[j][1];
				tnv.iz = tz + nb[j][2];
				var index = tnv.ix * pWidth * pHeight + pHeight * tnv.iy
						+ tnv.iz;

				if (tnv.ix < pLength && tnv.ix > -1 && tnv.iy < pWidth
						&& tnv.iy > -1 && tnv.iz < pHeight && tnv.iz > -1
						&& (vpBits[index] & INOUT) && !(vpBits[index] & ISDONE)) {
					boundPoint.set(tnv.ix, tnv.iy, tz + nb[j][2], bp);

					dx = tnv.ix - bp.ix;
					dy = tnv.iy - bp.iy;
					dz = tnv.iz - bp.iz;
					square = dx * dx + dy * dy + dz * dz;
					vpDistance[index] = Math.sqrt(square);
					vpBits[index] |= ISDONE;
					vpBits[index] |= ISBOUND;

					outarray.push({
						ix : tnv.ix,
						iy : tnv.iy,
						iz : tnv.iz
					});
					positout++;
				} else if (tnv.ix < pLength && tnv.ix > -1 && tnv.iy < pWidth
						&& tnv.iy > -1 && tnv.iz < pHeight && tnv.iz > -1
						&& (vpBits[index] & INOUT) && (vpBits[index] & ISDONE)) {
					dx = tnv.ix - bp.ix;
					dy = tnv.iy - bp.iy;
					dz = tnv.iz - bp.iz;
					square = Math.sqrt(dx * dx + dy * dy + dz * dz);
					if (square < vpDistance[index]) {
						boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);
						vpDistance[index] = square;
						if (!(vpBits[index] & ISBOUND)) {
							vpBits[index] |= ISBOUND;
							outarray.push({
								ix : tnv.ix,
								iy : tnv.iy,
								iz : tnv.iz
							});
							positout++;
						}
					}
				}
			}
		}

		// console.log("part2", positout);

		for (i = 0; i < number; i++) {
			tx = inarray[i].ix;
			ty = inarray[i].iy;
			tz = inarray[i].iz;
			var bp = boundPoint.get(tx, ty, tz);

			for (j = 18; j < 26; j++) {
				tnv.ix = tx + nb[j][0];
				tnv.iy = ty + nb[j][1];
				tnv.iz = tz + nb[j][2];
				var index = tnv.ix * pWidth * pHeight + pHeight * tnv.iy
						+ tnv.iz;

				if (tnv.ix < pLength && tnv.ix > -1 && tnv.iy < pWidth
						&& tnv.iy > -1 && tnv.iz < pHeight && tnv.iz > -1
						&& (vpBits[index] & INOUT) && !(vpBits[index] & ISDONE)) {
					boundPoint.set(tnv.ix, tnv.iy, tz + nb[j][2], bp);

					dx = tnv.ix - bp.ix;
					dy = tnv.iy - bp.iy;
					dz = tnv.iz - bp.iz;
					square = dx * dx + dy * dy + dz * dz;
					vpDistance[index] = Math.sqrt(square);
					vpBits[index] |= ISDONE;
					vpBits[index] |= ISBOUND;

					outarray.push({
						ix : tnv.ix,
						iy : tnv.iy,
						iz : tnv.iz
					});
					positout++;
				} else if (tnv.ix < pLength && tnv.ix > -1 && tnv.iy < pWidth
						&& tnv.iy > -1 && tnv.iz < pHeight && tnv.iz > -1
						&& (vpBits[index] & INOUT) && (vpBits[index] & ISDONE)) {
					dx = tnv.ix - bp.ix;
					dy = tnv.iy - bp.iy;
					dz = tnv.iz - bp.iz;
					square = Math.sqrt(dx * dx + dy * dy + dz * dz);
					if (square < vpDistance[index]) {
						boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);

						vpDistance[index] = square;
						if (!(vpBits[index] & ISBOUND)) {
							vpBits[index] |= ISBOUND;
							outarray.push({
								ix : tnv.ix,
								iy : tnv.iy,
								iz : tnv.iz
							});
							positout++;
						}
					}
				}
			}
		}

		// console.log("part3", positout);
		return positout;
	};

	this.marchingcubeinit = function(stype) {
		for ( var i = 0, lim = vpBits.length; i < lim; i++) {
			if (stype == 3) {// vdw
				vpBits[i] &= ~ISBOUND;
			} else if (stype == 4) { // ses
				vpBits[i] &= ~ISDONE;
				if (vpBits[i] & ISBOUND)
					vpBits[i] |= ISDONE;
				vpBits[i] &= ~ISBOUND;
			} else if (stype == 2) {// after vdw
				if ((vpBits[i] & ISBOUND) && (vpBits[i] & ISDONE))
					vpBits[i] &= ~ISBOUND;
				else if ((vpBits[i] & ISBOUND) && !(vpBits[i] & ISDONE))
					vpBits[i] |= ISDONE;
			} else if (stype == 3) { // sas
				vpBits[i] &= ~ISBOUND;
			}
		}
	};

	// create (or retrieve) a vertex at the appropriate point for
	// the edge (p1,p2)
	var getVertex = function(i, j, k, code, p1, p2, vertnums) {
		var val1 = !!(code & (1 << p1));
		var val2 = !!(code & (1 << p2));

		// p1 if they are the same or if !val1
		var p = p1;
		if (!val1 && val2)
			p = p2;

		// adjust i,j,k by p
		if (p & 1)
			k++;
		if (p & 2)
			j++;
		if (p & 4)
			i++;

		var index = ((pWidth * i) + j) * pHeight + k;
		if (vertnums[index] < 0) // not created yet
		{
			vertnums[index] = verts.length;
			verts.push({
				x : i,
				y : j,
				z : k
			});
		}
		return vertnums[index];
	};

	// this code allows me to empirically prune the marching cubes code tables
	// to more efficiently handle discrete data
	var counter = function() {
		var data = Array(256);
		for ( var i = 0; i < 256; i++)
			data[i] = [];

		this.incrementUsed = function(i, j) {
			if (typeof data[i][j] === 'undefined')
				data[i][j] = {
					used : 0,
					unused : 0
				};
			data[i][j].used++;
		};

		this.incrementUnused = function(i, j) {
			if (typeof data[i][j] === 'undefined')
				data[i][j] = {
					used : 0,
					unused : 0
				};
			data[i][j].unused++;

		};

		var redoTable = function(triTable) {
			var str = "[";
			for ( var i = 0; i < triTable.length; i++) {
				var code = 0;
				var table = triTable[i];
				for ( var j = 0; j < table.length; j++) {
					code |= (1 << (table[j]));
				}
				str += "0x" + code.toString(16) + ", ";
			}
			str += "]";
			console.log(str);
		}

		this.print = function() {

			var table = MarchingCube.triTable;
			var str;
			var newtable = [];
			for ( var i = 0; i < table.length; i++) {
				var newarr = [];
				for ( var j = 0; j < table[i].length; j += 3) {
					var k = j / 3;
					if (typeof data[i][k] === 'undefined' || !data[i][k].unused) {
						newarr.push(table[i][j]);
						newarr.push(table[i][j + 1]);
						newarr.push(table[i][j + 2]);
					}
					if (typeof data[i][k] === 'undefined')
						console.log("undef " + i + "," + k);
				}
				newtable.push(newarr);
			}
			console.log(JSON.stringify(newtable));
			redoTable(newtable);
		};
	};
	// this is based off the code here:
	// http://paulbourke.net/geometry/polygonise/
	// which is in turn based off of assorted public domain codes
	this.marchingcube = function(stype) {
		this.marchingcubeinit(stype);

		// this array keeps track of unique numbers for vertices
		// created in verts
		var vertnums = new Int32Array(pLength * pWidth * pHeight);
		for ( var i = 0; i < vertnums.length; i++) {
			vertnums[i] = -1;
		}

		verts = new Array();
		faces = new Array();

		// consider every grid cube
		var etable = MarchingCube.edgeTable;
		var tritable = MarchingCube.triTable;
		var i, j, k, p, t;
		var l, w, h, trilen, vlen;
		var vertList = new Int32Array(12);
		for (i = 0, l = pLength - 1; i < l; i++) {
			for (j = 0, w = pWidth - 1; j < w; j++) {
				for (k = 0, h = pHeight - 1; k < h; k++) {
					var code = 0;
					for (p = 0; p < 8; p++) {
						var index = ((pWidth * (i + ((p & 4) >> 2))) + j + ((p & 2) >> 1))
								* pHeight + k + (p & 1);

						var val = !!(vpBits[index] & ISDONE);
						code |= val << p;
					}

					// set the vertList
					var ecode = etable[code];
					if (ecode == 0)
						continue;
					var ttable = tritable[code];

					// based on code, determine what vertices are needed for
					// cube i,j,k
					// if necessary create them (adding to verts and setting
					// vertnums)
					// and set the appropriate vertex indices in vertList
					if (ecode & 1)
						vertList[0] = getVertex(i, j, k, code, 0, 1, vertnums);
					if (ecode & 2)
						vertList[1] = getVertex(i, j, k, code, 1, 3, vertnums);
					if (ecode & 4)
						vertList[2] = getVertex(i, j, k, code, 3, 2, vertnums);
					if (ecode & 8)
						vertList[3] = getVertex(i, j, k, code, 2, 0, vertnums);
					if (ecode & 16)
						vertList[4] = getVertex(i, j, k, code, 4, 5, vertnums);
					if (ecode & 32)
						vertList[5] = getVertex(i, j, k, code, 5, 7, vertnums);
					if (ecode & 64)
						vertList[6] = getVertex(i, j, k, code, 7, 6, vertnums);
					if (ecode & 128)
						vertList[7] = getVertex(i, j, k, code, 6, 4, vertnums);
					if (ecode & 256)
						vertList[8] = getVertex(i, j, k, code, 0, 4, vertnums);
					if (ecode & 512)
						vertList[9] = getVertex(i, j, k, code, 1, 5, vertnums);
					if (ecode & 1024)
						vertList[10] = getVertex(i, j, k, code, 3, 7, vertnums);
					if (ecode & 2048)
						vertList[11] = getVertex(i, j, k, code, 2, 6, vertnums);

					// add all faces
					for (t = 0, trilen = ttable.length; t < trilen; t += 3) {
						var a = vertList[ttable[t]];
						var b = vertList[ttable[t + 1]];
						var c = vertList[ttable[t + 2]];

						if (a != b && b != c && a != c)
							faces.push({
								a : a,
								b : b,
								c : c
							});
						/*
						 * if (a != b && b != c && a != c) {
						 * counts.incrementUsed(code, t / 3); } else
						 * counts.incrementUnused(code, t / 3);
						 */
					}

				}
			}
		}
		// set atom ids
		for (i = 0, vlen = verts.length; i < vlen; i++) {
			verts[i].atomid = vpAtomID[verts[i].x * pWidth * pHeight + pHeight
					* verts[i].y + verts[i].z];
		}
	};

});
self.onmessage = function(oEvent) {
	var obj = oEvent.data;
	var type = obj.type;

	if (type < 0) // sending atom data, initialize
	{
		self.atomData = obj.atoms;
		self.ps = new ProteinSurface();
	} else {
		var ps = self.ps;
		ps.initparm(obj.expandedExtent, (type == 1) ? false : true);
		ps.fillvoxels(self.atomData, obj.extendedAtoms);
		ps.buildboundary();

		if (type == 4 || type == 2)
			ps.fastdistancemap();
		if (type == 2) {
			ps.boundingatom(false);
			ps.fillvoxelswaals(self.atomData, obj.extendedAtoms);
		}

		ps.marchingcube(type);
		ps.laplaciansmooth(1);
		var VandF = ps.getFacesAndVertices(obj.atomsToShow);
		self.postMessage(VandF);
	}
};
