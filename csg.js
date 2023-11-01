// Ported from https://github.com/evanw/csg.js
import * as gl from './glConsts';
import { v3Negate, v3Dot, v3Cross, v3Sub, v3Normalize, vecLerp, v3Max, v3Length, v3Abs, v3AddScale, m4RotX, m4RotY, m4RotZ, m4Mul, m4MulPoint, v3Add, v3Mul } from "./types";
const CSG_PLANE_EPSILON = 1e-5;
let v3Scratch;
let csgPolygonNew = (verts, tag) => (v3Scratch = v3Normalize(v3Cross(v3Sub(verts[1].pos, verts[0].pos), v3Sub(verts[2].pos, verts[0].pos))),
    {
        vertices: verts,
        plane: { normal: v3Scratch, w: v3Dot(v3Scratch, verts[0].pos) },
        tag
    });
let csgPlaneFlip = (me) => ({
    normal: v3Negate(me.normal),
    w: -me.w,
});
let csgPlaneSplitPolygon = (me, polygon, coplanarFront, coplanarBackk, front, backk) => {
    let polygonType = 0;
    let types = polygon.vertices.map(vert => {
        let t = v3Dot(me.normal, vert.pos) - me.w;
        let typ = t < -CSG_PLANE_EPSILON ? 2 /* BACKK */
            : t > CSG_PLANE_EPSILON ? 1 /* FRONT */
                : 0 /* COPLANAR */;
        polygonType |= typ;
        return typ;
    });
    if (polygonType == 0 /* COPLANAR */) {
        (v3Dot(me.normal, polygon.plane.normal) > 0 ? coplanarFront : coplanarBackk).push(polygon);
    }
    if (polygonType == 1 /* FRONT */) {
        front.push(polygon);
    }
    if (polygonType == 2 /* BACKK */) {
        backk.push(polygon);
    }
    if (polygonType == 3 /* SPANNING */) {
        let f = [], b = [];
        for (let i = 0; i < polygon.vertices.length; i++) {
            let j = (i + 1) % polygon.vertices.length;
            let ti = types[i], tj = types[j];
            let vi = polygon.vertices[i], vj = polygon.vertices[j];
            if (ti != 2 /* BACKK */)
                f.push(vi);
            if (ti != 1 /* FRONT */)
                b.push(vi);
            if ((ti | tj) == 3 /* SPANNING */) {
                let t = (me.w - v3Dot(me.normal, vi.pos)) / v3Dot(me.normal, v3Sub(vj.pos, vi.pos));
                let v = {
                    pos: vecLerp(vi.pos, vj.pos, t),
                    normal: v3Normalize(vecLerp(vi.normal, vj.normal, t)),
                };
                f.push(v);
                b.push(v);
            }
        }
        if (f.length >= 3)
            front.push(csgPolygonNew(f, polygon.tag));
        if (b.length >= 3)
            backk.push(csgPolygonNew(b, polygon.tag));
    }
};
let csgNodeNew = () => ({
    plane: 0,
    front: 0,
    backk: 0,
    polygons: [],
});
let csgNodeInvert = (me) => {
    me.polygons = me.polygons.map(poly => ({
        vertices: poly.vertices.map(v => ({ pos: v.pos, normal: v3Negate(v.normal) })).reverse(),
        plane: csgPlaneFlip(poly.plane),
        tag: poly.tag,
    }));
    me.plane = csgPlaneFlip(me.plane); // assume not null
    if (me.front)
        csgNodeInvert(me.front);
    if (me.backk)
        csgNodeInvert(me.backk);
    let swap = me.front;
    me.front = me.backk;
    me.backk = swap;
};
let csgNodeClipPolygons = (me, polygons) => {
    if (!me.plane) {
        return [...me.polygons];
    }
    let front = [];
    let backk = [];
    polygons.map(poly => csgPlaneSplitPolygon(me.plane, poly, front, backk, front, backk));
    if (me.front) {
        front = csgNodeClipPolygons(me.front, front);
    }
    backk = me.backk
        ? csgNodeClipPolygons(me.backk, backk)
        : [];
    return [...front, ...backk];
};
let csgNodeClipTo = (me, bsp) => {
    me.polygons = csgNodeClipPolygons(bsp, me.polygons);
    if (me.front)
        csgNodeClipTo(me.front, bsp);
    if (me.backk)
        csgNodeClipTo(me.backk, bsp);
};
let csgNodeAllPolygons = (me) => [
    ...me.polygons,
    ...(me.front ? csgNodeAllPolygons(me.front) : []),
    ...(me.backk ? csgNodeAllPolygons(me.backk) : []),
];
let csgNodeBuild = (me, polygons) => {
    if (polygons.length) {
        if (!me.plane) {
            me.plane = polygons[0].plane;
        }
        let front = [];
        let backk = [];
        polygons.map(poly => csgPlaneSplitPolygon(me.plane, poly, me.polygons, me.polygons, front, backk));
        if (front.length) {
            if (!me.front) {
                me.front = csgNodeNew();
            }
            csgNodeBuild(me.front, front);
        }
        if (backk.length) {
            if (!me.backk) {
                me.backk = csgNodeNew();
            }
            csgNodeBuild(me.backk, backk);
        }
    }
};
// ----------------------------------------------------------------------------
const V_POSITION = '__a';
const F_BOX = '__b';
const F_LINE = '__c';
const F_MIN = '__d';
const F_MAX = '__e';
const F_NEG = '__f';
export let csgSolidOpUnion = (solidA, solidB) => {
    let a = csgNodeNew(), b = csgNodeNew();
    csgNodeBuild(a, solidA.polys);
    csgNodeBuild(b, solidB.polys);
    csgNodeClipTo(a, b);
    csgNodeClipTo(b, a);
    csgNodeInvert(b);
    csgNodeClipTo(b, a);
    csgNodeInvert(b);
    csgNodeBuild(a, csgNodeAllPolygons(b));
    return EDITOR ? {
        polys: csgNodeAllPolygons(a),
        lineViewPolys: (solidA.lineViewPolys || solidA.polys.map(x => (x.tag = 0, x)))
            .concat(solidB.lineViewPolys || solidB.polys.map(x => (x.tag = 0, x))),
        sdf: `${F_MIN}(${solidA.sdf},${solidB.sdf})`,
    } : {
        polys: csgNodeAllPolygons(a),
        sdf: `${F_MIN}(${solidA.sdf},${solidB.sdf})`,
    };
};
export let csgSolidOpSubtract = (solidA, solidB) => {
    let a = csgNodeNew(), b = csgNodeNew();
    csgNodeBuild(a, solidA.polys);
    csgNodeBuild(b, solidB.polys);
    csgNodeInvert(a);
    csgNodeClipTo(a, b);
    csgNodeClipTo(b, a);
    csgNodeInvert(b);
    csgNodeClipTo(b, a);
    csgNodeInvert(b);
    csgNodeBuild(a, csgNodeAllPolygons(b));
    csgNodeInvert(a);
    return EDITOR ? {
        polys: csgNodeAllPolygons(a),
        lineViewPolys: (solidA.lineViewPolys || solidA.polys.map(x => (x.tag = 0, x)))
            .concat((solidB.lineViewPolys || solidB.polys).map(x => (x.tag = 1, x))),
        sdf: `${F_MAX}(${solidA.sdf},${F_NEG}(${solidB.sdf}))`,
    } : {
        polys: csgNodeAllPolygons(a),
        sdf: `${F_MAX}(${solidA.sdf},${F_NEG}(${solidB.sdf}))`,
    };
};
export let csgSolidOpIntersect = (solidA, solidB) => {
    let a = csgNodeNew(), b = csgNodeNew();
    csgNodeBuild(a, solidA.polys);
    csgNodeBuild(b, solidB.polys);
    csgNodeInvert(a);
    csgNodeClipTo(b, a);
    csgNodeInvert(b);
    csgNodeClipTo(a, b);
    csgNodeClipTo(b, a);
    csgNodeBuild(a, csgNodeAllPolygons(b));
    csgNodeInvert(a);
    return EDITOR ? {
        polys: csgNodeAllPolygons(a),
        lineViewPolys: (solidA.lineViewPolys || solidA.polys.map(x => (x.tag = 0, x)))
            .concat((solidB.lineViewPolys || solidB.polys).map(x => (x.tag = 1, x))),
        sdf: ''
    } : {
        polys: csgNodeAllPolygons(a),
        sdf: ''
    };
};
let rot;
let accVertices;
let sphereVertexCenter;
let sphereVertexOffsetScale;
// theta: 0-4 longitude, phi: 0-2 latitude pole to pole
let sphereVertex = (radius, offset, theta, phi) => {
    theta *= Math.PI / 2;
    phi *= Math.PI / 2;
    let normal = [
        Math.cos(theta) * Math.sin(phi),
        Math.cos(phi),
        Math.sin(theta) * Math.sin(phi),
    ];
    accVertices.push({
        pos: v3Add(sphereVertexCenter, m4MulPoint(rot, v3AddScale(v3Mul(offset, sphereVertexOffsetScale), normal, radius))),
        normal: m4MulPoint(rot, normal),
    });
};
export let csgSolidLine = (tag, cx, cy, cz, h, r0, r1, yaw, pitch, roll) => {
    //let maxRad = Math.max(r0, r1)
    let resolution = 4; //Math.max(4,Math.round(maxRad/200))
    let polys = [];
    sphereVertexCenter = [cx, cy, cz];
    sphereVertexOffsetScale = [1, 1, 1];
    rot = m4Mul(m4Mul(m4RotY(yaw / 180 * Math.PI), m4RotX(pitch / 180 * Math.PI)), m4RotZ(roll / 180 * Math.PI));
    let d = r0 - r1;
    let phi = Math.atan(d / h);
    let phiLat = 2 * phi / Math.PI;
    let offsetY0 = r0 * Math.sin(phi);
    let walkR0 = r0 * Math.cos(phi);
    let offsetY1 = h + r1 * Math.sin(phi);
    let walkR1 = r1 * Math.cos(phi);
    for (let i = 0; i < 4 * resolution; ++i) { // longitudes
        let i0 = i / resolution, i1 = (i + 1) / resolution;
        accVertices = [];
        sphereVertex(walkR0, [0, offsetY0, 0], i0, 1);
        sphereVertex(walkR1, [0, offsetY1, 0], i0, 1);
        sphereVertex(walkR1, [0, offsetY1, 0], i1, 1);
        sphereVertex(walkR0, [0, offsetY0, 0], i1, 1);
        polys.push(csgPolygonNew(accVertices, tag));
        // top latitudes
        for (let j = 0; j < 2 * resolution; ++j) {
            let j0 = j / resolution, j1 = (j + 1) / resolution;
            let brk = j1 > 1 - phiLat;
            if (brk)
                j1 = 1 - phiLat;
            accVertices = [];
            sphereVertex(r1, [0, h, 0], i0, j0);
            j0 > 0.01 && sphereVertex(r1, [0, h, 0], i1, j0);
            sphereVertex(r1, [0, h, 0], i1, j1);
            sphereVertex(r1, [0, h, 0], i0, j1);
            polys.push(csgPolygonNew(accVertices, tag));
            if (brk)
                break;
        }
        // bottom latitudes
        for (let j = 2 * resolution; j > 0; --j) {
            let j0 = (j - 1) / resolution, j1 = j / resolution;
            let brk = j0 < 1 - phiLat;
            if (brk)
                j0 = 1 - phiLat;
            accVertices = [];
            sphereVertex(r0, [0, 0, 0], i0, j0);
            sphereVertex(r0, [0, 0, 0], i1, j0);
            j1 < 1.99 && sphereVertex(r0, [0, 0, 0], i1, j1);
            sphereVertex(r0, [0, 0, 0], i0, j1);
            polys.push(csgPolygonNew(accVertices, tag));
            if (brk)
                break;
        }
    }
    return {
        polys,
        sdf: `${F_LINE}(${tag},${V_POSITION},[${cx},${cy},${cz}],${h},${r0},${r1},${yaw},${pitch},${roll})`,
    };
};
export let csgSolidBox = (tag, cx, cy, cz, rx, ry, rz, yaw, pitch, roll, radius) => {
    let resolution = Math.max(4, Math.round(radius / 100));
    let polys = [];
    rot = m4Mul(m4Mul(m4RotY(yaw / 180 * Math.PI), m4RotX(pitch / 180 * Math.PI)), m4RotZ(roll / 180 * Math.PI));
    sphereVertexCenter = [cx, cy, cz];
    sphereVertexOffsetScale = [-rx, -ry, -rz];
    let defEdge = (i0, i1, offsets0, offsets1, mulA, addA, mulB, addB) => {
        let offset0 = [...offsets0].map(x => 2 * x - 1);
        let offset1 = [...offsets1].map(x => 2 * x - 1);
        accVertices = [];
        sphereVertex(radius, offset0, i0 * mulA + addA, i0 * mulB + addB);
        sphereVertex(radius, offset1, i0 * mulA + addA, i0 * mulB + addB);
        sphereVertex(radius, offset1, i1 * mulA + addA, i1 * mulB + addB);
        sphereVertex(radius, offset0, i1 * mulA + addA, i1 * mulB + addB);
        polys.push(csgPolygonNew(accVertices, tag));
    };
    if (radius > 0) {
        // corners
        for (let i = 0; i < resolution; ++i) { // longitudes
            for (let j = 0; j < resolution; ++j) { // latitudes
                for (let k = 0; k < 8; ++k) { // corner
                    let i0 = i / resolution + k % 4, i1 = (i + 1) / resolution + k % 4;
                    let j0 = j / resolution + (k / 4 | 0), j1 = (j + 1) / resolution + (k / 4 | 0);
                    let offset = [
                        (2 * (!!(k & 1) ^ !!(k & 2)) - 1),
                        (2 * !!(k & 4) - 1),
                        (2 * !!(k & 2) - 1),
                    ];
                    accVertices = [];
                    sphereVertex(radius, offset, i0, j0);
                    j0 > 0.01 && sphereVertex(radius, offset, i1, j0);
                    j1 < 1.99 && sphereVertex(radius, offset, i1, j1);
                    sphereVertex(radius, offset, i0, j1);
                    polys.push(csgPolygonNew(accVertices, tag));
                    // edges
                    if (!k && !j) {
                        defEdge(i0, i1, '011', '010', 0, 0, 1, 1);
                        defEdge(i0, i1, '110', '111', 0, 2, 1, 1);
                        defEdge(i0, i1, '001', '000', 0, 0, 1, 0);
                        defEdge(i0, i1, '100', '101', 0, 2, 1, 0);
                        defEdge(i0, i1, '010', '110', 0, 1, 1, 1);
                        defEdge(i0, i1, '111', '011', 0, 3, 1, 1);
                        defEdge(i0, i1, '000', '100', 0, 1, 1, 0);
                        defEdge(i0, i1, '101', '001', 0, 3, 1, 0);
                        defEdge(i0, i1, '010', '000', 1, 0, 0, 1);
                        defEdge(i0, i1, '110', '100', 1, 1, 0, 1);
                        defEdge(i0, i1, '111', '101', 1, 2, 0, 1);
                        defEdge(i0, i1, '011', '001', 1, 3, 0, 1);
                    }
                }
            }
        }
    }
    if (rx && ry && rz) {
        // faces
        [
            [[0, 4, 6, 2], [-1, 0, 0]],
            [[1, 3, 7, 5], [1, 0, 0]],
            [[0, 1, 5, 4], [0, -1, 0]],
            [[2, 6, 7, 3], [0, 1, 0]],
            [[0, 2, 3, 1], [0, 0, -1]],
            [[4, 5, 7, 6], [0, 0, 1]]
        ].map(info => polys.push(csgPolygonNew(info[0].map(i => (v3Scratch = v3AddScale([
            rx * (2 * !!(i & 1) - 1),
            ry * (2 * !!(i & 2) - 1),
            rz * (2 * !!(i & 4) - 1)
        ], info[1], radius),
            {
                pos: v3Add([cx, cy, cz], m4MulPoint(rot, v3Scratch)),
                normal: m4MulPoint(rot, info[1]),
            })), tag)));
    }
    return {
        polys,
        sdf: `${F_BOX}(${tag},${V_POSITION},[${cx},${cy},${cz}],[${rx},${ry},${rz}],${yaw},${pitch},${roll},${radius})`
    };
};
let sdfBox = (tag, p, center, extents, yaw, pitch, roll, radius) => (v3Scratch = v3Sub(v3Abs(m4MulPoint(m4Mul(m4Mul(m4RotZ(-roll / 180 * Math.PI), m4RotX(-pitch / 180 * Math.PI)), m4RotY(-yaw / 180 * Math.PI)), v3Sub(p, center))), extents),
    [tag, v3Length(v3Max(v3Scratch, [0, 0, 0])) + Math.min(Math.max(...v3Scratch), 0) - radius]);
let sdfLine = (tag, p, center, h, r0, r1, yaw, pitch, roll) => {
    v3Scratch = m4MulPoint(m4Mul(m4Mul(m4RotZ(-roll / 180 * Math.PI), m4RotX(-pitch / 180 * Math.PI)), m4RotY(-yaw / 180 * Math.PI)), v3Sub(p, center));
    let b = (r0 - r1) / h;
    let a = Math.sqrt(1.0 - b * b);
    let q = [v3Length([v3Scratch[0], v3Scratch[2], 0]), v3Scratch[1], 0];
    let k = v3Dot(q, [-b, a, 0]);
    return [tag, k < 0.0 ? v3Length(q) - r0 :
            k > a * h ? v3Length(v3Sub(q, [0, h, 0])) - r1 :
                v3Dot(q, [a, b, 0]) - r0];
};
let sdfMin = (a, b) => a[1] <= b[1] ? a : b;
let sdfMax = (a, b) => a[1] > b[1] ? a : b;
let sdfNeg = (a) => [a[0], -a[1]];
export let csgSolidBake = (me) => {
    let vertexBuf = [];
    let normalBuf = [];
    let tagBuf = [];
    let indexBuf = [];
    let linesIndexBuf = [];
    let linesVertexBuf = [];
    let linesTagBuf = [];
    let innerSdfFunc = new Function(`${V_POSITION},${F_BOX},${F_LINE},${F_MIN},${F_MAX},${F_NEG}`, 'return ' + me.sdf);
    let sdfFunc = (x) => innerSdfFunc(x, sdfBox, sdfLine, sdfMin, sdfMax, sdfNeg);
    me.polys.map(poly => {
        let startIdx = vertexBuf.length / 3;
        poly.vertices.map(x => (vertexBuf.push(...x.pos),
            normalBuf.push(...x.normal),
            tagBuf.push(poly.tag)));
        for (let i = 2; i < poly.vertices.length; i++) {
            indexBuf.push(startIdx, startIdx + i - 1, startIdx + i);
        }
    });
    if (EDITOR) {
        let showPolys = window.editorShowLinesKind ? me.polys : (me.lineViewPolys || me.polys);
        showPolys.map(poly => {
            let startIdx = linesVertexBuf.length / 3;
            poly.vertices.map(x => (linesVertexBuf.push(...x.pos),
                linesTagBuf.push(poly.tag)));
            for (let i = 2; i < poly.vertices.length; i++) {
                linesIndexBuf.push(startIdx, startIdx + i - 1);
                linesIndexBuf.push(startIdx + i - 1, startIdx + i);
                linesIndexBuf.push(startIdx + i, startIdx);
            }
        });
    }
    let index = G.createBuffer();
    G.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, index);
    G.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(indexBuf), gl.STATIC_DRAW);
    let vertex = G.createBuffer();
    G.bindBuffer(gl.ARRAY_BUFFER, vertex);
    G.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertexBuf), gl.STATIC_DRAW);
    let normal = G.createBuffer();
    G.bindBuffer(gl.ARRAY_BUFFER, normal);
    G.bufferData(gl.ARRAY_BUFFER, new Float32Array(normalBuf), gl.STATIC_DRAW);
    let tag = G.createBuffer();
    G.bindBuffer(gl.ARRAY_BUFFER, tag);
    G.bufferData(gl.ARRAY_BUFFER, new Float32Array(tagBuf), gl.STATIC_DRAW);
    let linesIndex;
    let linesVertex;
    let linesTag;
    if (EDITOR) {
        linesIndex = G.createBuffer();
        G.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, linesIndex);
        G.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(linesIndexBuf), gl.STATIC_DRAW);
        linesVertex = G.createBuffer();
        G.bindBuffer(gl.ARRAY_BUFFER, linesVertex);
        G.bufferData(gl.ARRAY_BUFFER, new Float32Array(linesVertexBuf), gl.STATIC_DRAW);
        linesTag = G.createBuffer();
        G.bindBuffer(gl.ARRAY_BUFFER, linesTag);
        G.bufferData(gl.ARRAY_BUFFER, new Float32Array(linesTagBuf), gl.STATIC_DRAW);
    }
    return [
        EDITOR ? {
            indexBuffer: index,
            indexBufferLen: indexBuf.length,
            vertexBuffer: vertex,
            normalBuffer: normal,
            tagBuffer: tag,
            lines: {
                indexBuffer: linesIndex,
                indexBufferLen: linesIndexBuf.length,
                vertexBuffer: linesVertex,
                tagBuffer: linesTag,
            }
        } : {
            indexBuffer: index,
            indexBufferLen: indexBuf.length,
            vertexBuffer: vertex,
            normalBuffer: normal,
            tagBuffer: tag,
        },
        sdfFunc
    ];
};
export let modelGeoDelete = (geo) => {
    if (!EDITOR)
        return;
    G.deleteBuffer(geo.indexBuffer);
    G.deleteBuffer(geo.vertexBuffer);
    G.deleteBuffer(geo.normalBuffer);
    G.deleteBuffer(geo.tagBuffer);
    if (geo.lines) {
        G.deleteBuffer(geo.lines.indexBuffer);
        G.deleteBuffer(geo.lines.vertexBuffer);
        G.deleteBuffer(geo.lines.tagBuffer);
    }
};
