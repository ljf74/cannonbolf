const LINES = 1;
const TRIANGLES = 4;
const SRC_ALPHA = 770;
const ONE_MINUS_SRC_ALPHA = 771;
const ARRAY_BUFFER = 34962;
const ELEMENT_ARRAY_BUFFER = 34963;
const STATIC_DRAW = 35044;
const DYNAMIC_DRAW = 35048;
const TEXTURE_2D = 3553;
const CULL_FACE = 2884;
const BLEND = 3042;
const DEPTH_TEST = 2929;
const UNSIGNED_BYTE = 5121;
const UNSIGNED_SHORT = 5123;
const FLOAT = 5126;
const RGBA = 6408;
const FRAGMENT_SHADER = 35632;
const VERTEX_SHADER = 35633;
const LEQUAL = 515;
const NEAREST = 9728;
const TEXTURE_MAG_FILTER = 10240;
const TEXTURE_MIN_FILTER = 10241;
const TEXTURE_WRAP_S = 10242;
const TEXTURE_WRAP_T = 10243;
const TEXTURE0 = 33984;
const REPEAT = 10497;
const CLAMP_TO_EDGE = 33071;

const False = 0;
const True = 1;
let lerp = (a, b, t) => a + t * (b - a);
let v3Negate = (a) => a.map(x => -x);
let v3Add = (a, b) => a.map((x, i) => x + b[i]);
let v3Sub = (a, b) => a.map((x, i) => x - b[i]);
let v3Mul = (a, b) => a.map((x, i) => x * b[i]);
let v3AddScale = (a, b, s) => a.map((x, i) => x + s * b[i]);
let v3Abs = (a) => a.map(Math.abs);
let v3Max = (a, b) => a.map((x, i) => Math.max(x, b[i]));
let vecLerp = (a, b, t) => a.map((x, i) => lerp(x, b[i], t));
let radLerp = (a, b, t) => {
    let delta = b - a;
    let lerp = delta > Math.PI
        ? delta - 2 * Math.PI
        : delta < -Math.PI
            ? delta + 2 * Math.PI
            : delta;
    let ret = a + lerp * t;
    if (ret < 0)
        ret += 2 * Math.PI;
    if (ret > 2 * Math.PI)
        ret -= 2 * Math.PI;
    return ret;
};
let v3Dot = ([x, y, z], [a, b, c]) => x * a + y * b + z * c;
let v3Dot2 = (a) => v3Dot(a, a);
let v3Cross = ([x, y, z], [a, b, c]) => [y * c - z * b, z * a - x * c, x * b - y * a];
let v3Length = (x) => Math.sqrt(v3Dot2(x));
let v3Normalize = (a) => v3AddScale([0, 0, 0], a, 1 / (v3Length(a) || 1));
let v3Reflect = (v, norm, normScale, tanScale) => {
    let sign = Math.sign(norm[2]) || 1;
    let a = -1 / (sign + norm[2]);
    let b = norm[0] * norm[1] * a;
    let t0 = [1 + sign * norm[0] * norm[0] * a, sign * b, -sign * norm[0]];
    let t1 = [b, sign + norm[1] * norm[1] * a, -norm[1]];
    let vn = -normScale * v3Dot(v, norm);
    let vt0 = tanScale * v3Dot(v, t0);
    let vt1 = tanScale * v3Dot(v, t1);
    let ret = [0, 0, 0];
    ret = v3AddScale(ret, norm, vn);
    ret = v3AddScale(ret, t0, vt0);
    ret = v3AddScale(ret, t1, vt1);
    return ret;
};
let m4Perspective = (invAspect, near, far) => [
    // FOV = PI / 2
    invAspect, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, (far + near) / (near - far), -1,
    0, 0, (2 * far * near) / (near - far), 0
];
let c, s;
let m4RotX = (rads) => (c = Math.cos(rads), s = Math.sin(rads),
    [
        1, 0, 0, 0,
        0, c, s, 0,
        0, -s, c, 0,
        0, 0, 0, 1
    ]);
let m4RotY = (rads) => (c = Math.cos(rads), s = Math.sin(rads),
    [
        c, 0, s, 0,
        0, 1, 0, 0,
        -s, 0, c, 0,
        0, 0, 0, 1
    ]);
let m4RotZ = (rads) => (c = Math.cos(rads), s = Math.sin(rads),
    [
        c, s, 0, 0,
        -s, c, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    ]);
let m4Scale = (s) => ([
    s, 0, 0, 0,
    0, s, 0, 0,
    0, 0, s, 0,
    0, 0, 0, 1
]);
let m4Translate = ([x, y, z]) => [
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    x, y, z, 1
];
let m4AxisAngle = ([x, y, z], angleRads) => {
    let s = Math.sin(angleRads);
    let c = Math.cos(angleRads);
    let t = 1 - c;
    return [
        x * x * t + c,
        y * x * t + z * s,
        z * x * t - y * s,
        0,
        x * y * t - z * s,
        y * y * t + c,
        z * y * t + x * s,
        0,
        x * z * t + y * s,
        y * z * t - x * s,
        z * z * t + c,
        0,
        0, 0, 0, 1
    ];
};
let m4Ident = m4Scale(1);
let m4Mul = (a, b) => {
    let result = [];
    for (let x = 0; x < 16; ++x) {
        let i = 4 * (x / 4 | 0), j = x % 4;
        result.push(b[i] * a[j] + b[i + 1] * a[j + 4] + b[i + 2] * a[j + 8] + b[i + 3] * a[j + 12]);
    }
    return result;
};
let m4MulPoint = (m, [x, y, z]) => {
    let s = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1;
    return [
        (m[0] * x + m[4] * y + m[8] * z + m[12]) / s,
        (m[1] * x + m[5] * y + m[9] * z + m[13]) / s,
        (m[2] * x + m[6] * y + m[10] * z + m[14]) / s
    ];
};

// Ported from https://github.com/evanw/csg.js
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
let csgSolidOpUnion = (solidA, solidB) => {
    let a = csgNodeNew(), b = csgNodeNew();
    csgNodeBuild(a, solidA.polys);
    csgNodeBuild(b, solidB.polys);
    csgNodeClipTo(a, b);
    csgNodeClipTo(b, a);
    csgNodeInvert(b);
    csgNodeClipTo(b, a);
    csgNodeInvert(b);
    csgNodeBuild(a, csgNodeAllPolygons(b));
    return {
        polys: csgNodeAllPolygons(a),
        sdf: `${F_MIN}(${solidA.sdf},${solidB.sdf})`,
    };
};
let csgSolidOpSubtract = (solidA, solidB) => {
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
    return {
        polys: csgNodeAllPolygons(a),
        sdf: `${F_MAX}(${solidA.sdf},${F_NEG}(${solidB.sdf}))`,
    };
};
let csgSolidOpIntersect = (solidA, solidB) => {
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
    return {
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
let csgSolidLine = (tag, cx, cy, cz, h, r0, r1, yaw, pitch, roll) => {
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
let csgSolidBox = (tag, cx, cy, cz, rx, ry, rz, yaw, pitch, roll, radius) => {
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
let csgSolidBake = (me) => {
    let vertexBuf = [];
    let normalBuf = [];
    let tagBuf = [];
    let indexBuf = [];
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
    let index = G.createBuffer();
    G.bindBuffer(ELEMENT_ARRAY_BUFFER, index);
    G.bufferData(ELEMENT_ARRAY_BUFFER, new Uint16Array(indexBuf), STATIC_DRAW);
    let vertex = G.createBuffer();
    G.bindBuffer(ARRAY_BUFFER, vertex);
    G.bufferData(ARRAY_BUFFER, new Float32Array(vertexBuf), STATIC_DRAW);
    let normal = G.createBuffer();
    G.bindBuffer(ARRAY_BUFFER, normal);
    G.bufferData(ARRAY_BUFFER, new Float32Array(normalBuf), STATIC_DRAW);
    let tag = G.createBuffer();
    G.bindBuffer(ARRAY_BUFFER, tag);
    G.bufferData(ARRAY_BUFFER, new Float32Array(tagBuf), STATIC_DRAW);
    return [
        {
            indexBuffer: index,
            indexBufferLen: indexBuf.length,
            vertexBuffer: vertex,
            normalBuffer: normal,
            tagBuffer: tag,
        },
        sdfFunc
    ];
};

/* -*- mode: javascript; tab-width: 4; indent-tabs-mode: nil; -*-
*
* Copyright (c) 2011-2013 Marcus Geelnard
*
* This software is provided 'as-is', without any express or implied
* warranty. In no event will the authors be held liable for any damages
* arising from the use of this software.
*
* Permission is granted to anyone to use this software for any purpose,
* including commercial applications, and to alter it and redistribute it
* freely, subject to the following restrictions:
*
* 1. The origin of this software must not be misrepresented; you must not
*    claim that you wrote the original software. If you use this software
*    in a product, an acknowledgment in the product documentation would be
*    appreciated but is not required.
*
* 2. Altered source versions must be plainly marked as such, and must not be
*    misrepresented as being the original software.
*
* 3. This notice may not be removed or altered from any source
*    distribution.
*
*/
// Some general notes and recommendations:
//  * This code uses modern ECMAScript features, such as ** instead of
//    Math.pow(). You may have to modify the code to make it work on older
//    browsers.
//  * If you're not using all the functionality (e.g. not all oscillator types,
//    or certain effects), you can reduce the size of the player routine even
//    further by deleting the code.
//
// This music an the example track from the SoundBox player (https://sb.bitsnbites.eu/):
//   Zebrain by esaruoho, from the 4k Matlab intro Zebrain
//
// This music has been exported by SoundBox. You can use it with
// http://sb.bitsnbites.eu/player-small.js in your own product.
// See http://sb.bitsnbites.eu/demo.html for an example of how to
// use it in a demo.
let doGen$1 = () => {
    // Song data
    var song = {
        songData: [
            {
                i: [
                    1,
                    15,
                    164,
                    0,
                    1,
                    127,
                    128,
                    15,
                    0,
                    36,
                    0,
                    28,
                    87,
                    0,
                    0,
                    0,
                    0,
                    64,
                    1,
                    1,
                    2,
                    38,
                    128,
                    1,
                    3,
                    255,
                    2,
                    212,
                    8 // FX_DELAY_TIME
                ],
                // Patterns
                p: [1, 1, 2, 2, 3, 1, 3, 1, 1, 1, 2, 2, 3, 1, 3, 1, 1, 1, 2, 2, 3, 1, 3, 1, 1],
                // Columns
                c: [
                    { n: [135, 147, 142, 149, 147, 140],
                        f: [] },
                    { n: [130, 147, 152, 149, 144, 154],
                        f: [] },
                    { n: [135, 147, 142, 149, 147, 140, , , , , , , 161, 156, 151, 135, , , , , , , , , 138, 150, 145, 157, 131, 143, 135, 147],
                        f: [] }
                ]
            },
            {
                i: [
                    0,
                    105,
                    128,
                    0,
                    0,
                    255,
                    128,
                    33,
                    0,
                    7,
                    7,
                    20,
                    39,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    2,
                    41,
                    0,
                    0,
                    8,
                    0,
                    6,
                    0,
                    0 // FX_DELAY_TIME
                ],
                // Patterns
                p: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                // Columns
                c: [
                    { n: [, , , 99, , , 99, , , , , , , , , 111, , , , 99, , , 99],
                        f: [] }
                ]
            },
            {
                i: [
                    1,
                    55,
                    128,
                    0,
                    1,
                    81,
                    128,
                    0,
                    0,
                    0,
                    52,
                    6,
                    0,
                    0,
                    0,
                    0,
                    0,
                    243,
                    10,
                    1,
                    2,
                    135,
                    0,
                    0,
                    3,
                    147,
                    6,
                    230,
                    4 // FX_DELAY_TIME
                ],
                // Patterns
                p: [, , , , , , 2, , , , 1, , , , 2, , 3, 4, 1, 5, 6, 7, 2, 8, 3, 4, 1, 5, 6, 7, 2, 8],
                // Columns
                c: [
                    { n: [137],
                        f: [] },
                    { n: [, , 135],
                        f: [] },
                    { n: [159],
                        f: [] },
                    { n: [157],
                        f: [] },
                    { n: [152],
                        f: [] },
                    { n: [154],
                        f: [] },
                    { n: [149],
                        f: [] },
                    { n: [147],
                        f: [] }
                ]
            },
            {
                i: [
                    1,
                    100,
                    128,
                    0,
                    1,
                    201,
                    128,
                    34,
                    0,
                    0,
                    5,
                    6,
                    58,
                    0,
                    0,
                    0,
                    0,
                    195,
                    6,
                    1,
                    2,
                    135,
                    0,
                    0,
                    9,
                    147,
                    6,
                    0,
                    1 // FX_DELAY_TIME
                ],
                // Patterns
                p: [, , , , , , , , 1, 1, 3, 3, 2, 3, 2, 3, 1, 1, 3, 3, 2, 3, 2, 3, 1, 1, 1, 1],
                // Columns
                c: [
                    { n: [111, , , 111, , , 111, , 106, , , , , , , , 104, , , 104, , , 104],
                        f: [] },
                    { n: [111, , , 111, , , 111, , 106, , , , , , , , 104, , , 104, , , 104, , 109, , , 109, , , 109],
                        f: [] },
                    { n: [111, , , 111, , , 111, , 106, , , , , , , , 104, , , 104, , , 104, , 109, , 109, , 109, , 109],
                        f: [] }
                ]
            },
            {
                i: [
                    0,
                    255,
                    116,
                    64,
                    0,
                    255,
                    116,
                    0,
                    64,
                    0,
                    4,
                    6,
                    35,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    2,
                    14,
                    0,
                    19,
                    32,
                    0,
                    0,
                    93,
                    1 // FX_DELAY_TIME
                ],
                // Patterns
                p: [, , , , , , , , 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3],
                // Columns
                c: [
                    { n: [, , 147, , , , , , , , , , 147, , , , 147, , , , , , , , 147, , , , , , , , , , , , , , , , 123, , , 123],
                        f: [] },
                    { n: [147, , , 147, , , , , , , , , 147, , , , 147, , , , , , , , 147, , , , , , , , , , , , , , , , 123, , , 123],
                        f: [] },
                    { n: [147],
                        f: [] }
                ]
            },
            {
                i: [
                    1,
                    11,
                    128,
                    0,
                    1,
                    111,
                    128,
                    25,
                    0,
                    0,
                    72,
                    6,
                    115,
                    0,
                    0,
                    0,
                    0,
                    53,
                    6,
                    1,
                    2,
                    135,
                    0,
                    0,
                    3,
                    255,
                    6,
                    168,
                    9 // FX_DELAY_TIME
                ],
                // Patterns
                p: [, , , , , , , , 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8],
                // Columns
                c: [
                    { n: [147, 138, 145, 150, 159, 149, 142, 161],
                        f: [] },
                    { n: [, , , , , , , , 140, , 159, 147, , 149, , 154, , , , 161],
                        f: [] },
                    { n: [, , , , , , , 145, , 150, , , 149, , , 147, , 143, , , , 154, , , 142],
                        f: [] },
                    { n: [, , , , , , 145, , , 149, , , 142, , , 138, , , , , 137, , , 135, , , 142],
                        f: [] },
                    { n: [, , , , , , , , , , , 150, 142, 145, , , , , , , , , , , , , , 152, 145],
                        f: [] },
                    { n: [, , , , , , , , , , , , 135, 150, 142, 159],
                        f: [] },
                    { n: [, , , , 149, , 154, 145, , , , , , , , , , , , , 138, , 143],
                        f: [] },
                    { n: [, , , , 138, 147, , , , 143, , , 149, , 142, , 150, 154, 145, , , , 154, 159, 149],
                        f: [] }
                ]
            },
            {
                i: [
                    3,
                    0,
                    128,
                    0,
                    3,
                    68,
                    128,
                    0,
                    64,
                    218,
                    4,
                    4,
                    40,
                    0,
                    0,
                    0,
                    1,
                    55,
                    4,
                    1,
                    2,
                    67,
                    115,
                    124,
                    190,
                    67,
                    6,
                    39,
                    1 // FX_DELAY_TIME
                ],
                // Patterns
                p: [, , , , , , , , , , , , 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                // Columns
                c: [
                    { n: [, , , 140, , , , , , , , , 147, , , , , 140, , 140, , , , , , , , , 147],
                        f: [] }
                ]
            },
        ],
        rowLen: 6615,
        patternLen: 32,
        endPattern: 31,
        numChannels: 7 // Number of channels
    };
    var CPlayer = function () {
        //--------------------------------------------------------------------------
        // Private methods
        //--------------------------------------------------------------------------
        // Oscillators
        var osc_sin = function (value) {
            return Math.sin(value * 6.283184);
        };
        var osc_saw = function (value) {
            return 2 * (value % 1) - 1;
        };
        var osc_square = function (value) {
            return (value % 1) < 0.5 ? 1 : -1;
        };
        var osc_tri = function (value) {
            var v2 = (value % 1) * 4;
            if (v2 < 2)
                return v2 - 1;
            return 3 - v2;
        };
        var getnotefreq = function (n) {
            // 174.61.. / 44100 = 0.003959503758 (F3)
            return 0.003959503758 * (2 ** ((n - 128) / 12));
        };
        var createNote = function (instr, n, rowLen) {
            var osc1 = mOscillators[instr.i[0]], o1vol = instr.i[1], o1xenv = instr.i[3] / 32, osc2 = mOscillators[instr.i[4]], o2vol = instr.i[5], o2xenv = instr.i[8] / 32, noiseVol = instr.i[9], attack = instr.i[10] * instr.i[10] * 4, sustain = instr.i[11] * instr.i[11] * 4, release = instr.i[12] * instr.i[12] * 4, releaseInv = 1 / release, expDecay = -instr.i[13] / 16, arp = instr.i[14], arpInterval = rowLen * (2 ** (2 - instr.i[15]));
            var noteBuf = new Int32Array(attack + sustain + release);
            // Re-trig oscillators
            var c1 = 0, c2 = 0;
            // Local variables.
            var j, j2, e, rsample, o1t, o2t;
            // Generate one note (attack + sustain + release)
            for (j = 0, j2 = 0; j < attack + sustain + release; j++, j2++) {
                if (j2 >= 0) {
                    // Switch arpeggio note.
                    arp = (arp >> 8) | ((arp & 255) << 4);
                    j2 -= arpInterval;
                    // Calculate note frequencies for the oscillators
                    o1t = getnotefreq(n + (arp & 15) + instr.i[2] - 128);
                    o2t = getnotefreq(n + (arp & 15) + instr.i[6] - 128) * (1 + 0.0008 * instr.i[7]);
                }
                // Envelope
                e = 1;
                if (j < attack) {
                    e = j / attack;
                }
                else if (j >= attack + sustain) {
                    e = (j - attack - sustain) * releaseInv;
                    e = (1 - e) * (3 ** (expDecay * e));
                }
                // Oscillator 1
                c1 += o1t * e ** o1xenv;
                rsample = osc1(c1) * o1vol;
                // Oscillator 2
                c2 += o2t * e ** o2xenv;
                rsample += osc2(c2) * o2vol;
                // Noise oscillator
                if (noiseVol) {
                    rsample += (2 * Math.random() - 1) * noiseVol;
                }
                // Add to (mono) channel buffer
                noteBuf[j] = (80 * rsample * e) | 0;
            }
            return noteBuf;
        };
        //--------------------------------------------------------------------------
        // Private members
        //--------------------------------------------------------------------------
        // Array of oscillator functions
        var mOscillators = [
            osc_sin,
            osc_square,
            osc_saw,
            osc_tri
        ];
        // Private variables set up by init()
        var mSong, mLastRow, mCurrentCol, mNumWords, mMixBuf;
        //--------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------
        this.init = function (song) {
            // Define the song
            mSong = song;
            // Init iteration state variables
            mLastRow = song.endPattern;
            mCurrentCol = 0;
            // Prepare song info
            mNumWords = song.rowLen * song.patternLen * (mLastRow + 1) * 2;
            // Create work buffer (initially cleared)
            mMixBuf = new Int32Array(mNumWords);
        };
        //--------------------------------------------------------------------------
        // Public methods
        //--------------------------------------------------------------------------
        // Generate audio data for a single track
        this.generate = function () {
            // Local variables
            var i, j, p, row, col, n, cp, k, t, rsample, rowStartSample, f;
            // Put performance critical items in local variables
            var chnBuf = new Int32Array(mNumWords), instr = mSong.songData[mCurrentCol], rowLen = mSong.rowLen, patternLen = mSong.patternLen;
            // Clear effect state
            var low = 0, band = 0, high;
            var lsample, filterActive = false;
            // Clear note cache.
            var noteCache = [];
            // Patterns
            for (p = 0; p <= mLastRow; ++p) {
                cp = instr.p[p];
                // Pattern rows
                for (row = 0; row < patternLen; ++row) {
                    // Execute effect command.
                    var cmdNo = cp ? instr.c[cp - 1].f[row] : 0;
                    if (cmdNo) {
                        instr.i[cmdNo - 1] = instr.c[cp - 1].f[row + patternLen] || 0;
                        // Clear the note cache since the instrument has changed.
                        if (cmdNo < 17) {
                            noteCache = [];
                        }
                    }
                    // Put performance critical instrument properties in local variables
                    var oscLFO = mOscillators[instr.i[16]], lfoAmt = instr.i[17] / 512, lfoFreq = (2 ** (instr.i[18] - 9)) / rowLen, fxLFO = instr.i[19], fxFilter = instr.i[20], fxFreq = instr.i[21] * 43.23529 * 3.141592 / 44100, q = 1 - instr.i[22] / 255, dist = instr.i[23] * 1e-5, drive = instr.i[24] / 32, panAmt = instr.i[25] / 512, panFreq = 6.283184 * (2 ** (instr.i[26] - 9)) / rowLen, dlyAmt = instr.i[27] / 255, dly = instr.i[28] * rowLen & ~1; // Must be an even number
                    // Calculate start sample number for this row in the pattern
                    rowStartSample = (p * patternLen + row) * rowLen;
                    // Generate notes for this pattern row
                    for (col = 0; col < 4; ++col) {
                        n = cp ? instr.c[cp - 1].n[row + col * patternLen] : 0;
                        if (n) {
                            if (!noteCache[n]) {
                                noteCache[n] = createNote(instr, n, rowLen);
                            }
                            // Copy note from the note cache
                            var noteBuf = noteCache[n];
                            for (j = 0, i = rowStartSample * 2; j < noteBuf.length; j++, i += 2) {
                                chnBuf[i] += noteBuf[j];
                            }
                        }
                    }
                    // Perform effects for this pattern row
                    for (j = 0; j < rowLen; j++) {
                        // Dry mono-sample
                        k = (rowStartSample + j) * 2;
                        rsample = chnBuf[k];
                        // We only do effects if we have some sound input
                        if (rsample || filterActive) {
                            // State variable filter
                            f = fxFreq;
                            if (fxLFO) {
                                f *= oscLFO(lfoFreq * k) * lfoAmt + 0.5;
                            }
                            f = 1.5 * Math.sin(f);
                            low += f * band;
                            high = q * (rsample - band) - low;
                            band += f * high;
                            rsample = fxFilter == 3 ? band : fxFilter == 1 ? high : low;
                            // Distortion
                            if (dist) {
                                rsample *= dist;
                                rsample = rsample < 1 ? rsample > -1 ? osc_sin(rsample * .25) : -1 : 1;
                                rsample /= dist;
                            }
                            // Drive
                            rsample *= drive;
                            // Is the filter active (i.e. still audiable)?
                            filterActive = rsample * rsample > 1e-5;
                            // Panning
                            t = Math.sin(panFreq * k) * panAmt + 0.5;
                            lsample = rsample * (1 - t);
                            rsample *= t;
                        }
                        else {
                            lsample = 0;
                        }
                        // Delay is always done, since it does not need sound input
                        if (k >= dly) {
                            // Left channel = left + right[-p] * t
                            lsample += chnBuf[k - dly + 1] * dlyAmt;
                            // Right channel = right + left[-p] * t
                            rsample += chnBuf[k - dly] * dlyAmt;
                        }
                        // Store in stereo channel buffer (needed for the delay effect)
                        chnBuf[k] = lsample | 0;
                        chnBuf[k + 1] = rsample | 0;
                        // ...and add to stereo mix buffer
                        mMixBuf[k] += lsample | 0;
                        mMixBuf[k + 1] += rsample | 0;
                    }
                }
            }
            // Next iteration. Return progress (1.0 == done!).
            mCurrentCol++;
            return mCurrentCol / mSong.numChannels;
        };
        // Create a WAVE formatted Uint8Array from the generated audio data
        this.createWave = function () {
            // Create WAVE header
            var headerLen = 44;
            var l1 = headerLen + mNumWords * 2 - 8;
            var l2 = l1 - 36;
            var wave = new Uint8Array(headerLen + mNumWords * 2);
            wave.set([82, 73, 70, 70,
                l1 & 255, (l1 >> 8) & 255, (l1 >> 16) & 255, (l1 >> 24) & 255,
                87, 65, 86, 69, 102, 109, 116, 32, 16, 0, 0, 0, 1, 0, 2, 0,
                68, 172, 0, 0, 16, 177, 2, 0, 4, 0, 16, 0, 100, 97, 116, 97,
                l2 & 255, (l2 >> 8) & 255, (l2 >> 16) & 255, (l2 >> 24) & 255]);
            // Append actual wave data
            for (var i = 0, idx = headerLen; i < mNumWords; ++i) {
                // Note: We clamp here
                var y = mMixBuf[i];
                y = y < -32767 ? -32767 : (y > 32767 ? 32767 : y);
                wave[idx++] = y & 255;
                wave[idx++] = (y >> 8) & 255;
            }
            // Return the WAVE formatted typed array
            return wave;
        };
        // Get n samples of wave data at time t [s]. Wave data in range [-2,2].
        this.getData = function (t, n) {
            var i = 2 * Math.floor(t * 44100);
            var d = new Array(n);
            for (var j = 0; j < 2 * n; j += 1) {
                var k = i + j;
                d[j] = t > 0 && k < mMixBuf.length ? mMixBuf[k] / 32768 : 0;
            }
            return d;
        };
    };
    var player = new CPlayer();
    player.init(song);
    var interval = setInterval(function () {
        if (player.generate() >= 1) {
            clearInterval(interval);
            var wave = player.createWave();
            self.postMessage(wave, wave.buffer);
        }
    }, 10);
};
let code = `let a=(${doGen$1.toString()})();self.postMessage(a,[a.buffer])`;
let worker = new Worker(URL.createObjectURL(new Blob([code])));
let prom = new Promise(resolve => {
    worker.onmessage = e => resolve(e.data);
});
let started = false;
let playMusic = () => {
    if (started)
        return;
    started = true;
    prom.then(wave => {
        var audio = document.createElement("audio");
        audio.src = URL.createObjectURL(new Blob([wave], { type: "audio/wav" }));
        audio.volume = 0.8;
        audio.loop = true;
        audio.play();
    });
};

let inputsNew = () => (({
    mouseAccX: 0,
    mouseAccY: 0,
    keysDown: {},
}));
let frame$1 = inputsNew();
let clickedIn = {};
let lastMouseDx = 0;
let lastMouseDy = 0;
document.onmousemove = (e) => {
    let dx = e.movementX, dy = e.movementY;
    if ((dx * dx < 10000 || dx * lastMouseDx < 0) && (dy * dy < 10000 || dy * lastMouseDy < 0)) {
        frame$1.mouseAccX += lastMouseDx = dx;
        frame$1.mouseAccY += lastMouseDy = dy;
    }
};
document.onmousedown = (e) => {
    if (document.pointerLockElement !== CC) {
        clickedIn.a = 1;
        playMusic();
        CC.requestPointerLock();
    }
    else {
        frame$1.keysDown[e.button] = True;
    }
};
document.onmouseup = (e) => {
    frame$1.keysDown[e.button] = False;
};
document.onkeydown = (e) => {
    frame$1.keysDown[e.code[3]] = True;
    return false;
};
document.onkeyup = (e) => {
    frame$1.keysDown[e.code[3]] = False;
};
let inputsConsumeFrame = () => {
    let outFrame = frame$1;
    frame$1 = inputsNew();
    inputsAdd(frame$1, outFrame);
    frame$1.mouseAccX = frame$1.mouseAccY = 0;
    return document.pointerLockElement == CC ? outFrame : inputsNew();
};
let inputsAdd = (me, other) => {
    for (let k in other.keysDown) {
        me.keysDown[k] |= other.keysDown[k];
    }
    me.mouseAccX += other.mouseAccX;
    me.mouseAccY += other.mouseAccY;
};

/* File generated with Shader Minifier 1.2
 * http://www.ctrl-alt-test.fr
 */
let aimRay_frag = `precision highp float;uniform vec3 a0;void main(){gl_FragColor=vec4(a0,.25);}`;
let aimRay_vert = `attribute vec3 a1;uniform mat4 a2,a3;void main(){gl_Position=a2*a3*vec4(a1,1);}`;
let blit_frag = `precision highp float;uniform sampler2D a4;uniform vec4 a5;void main(){vec2 a=gl_FragCoord.xy/a5.xy;vec4 i=texture2D(a4,vec2(a.x,1.-a.y));gl_FragColor=mix(vec4(mix(vec3(0),i.xyz,a5.w),1),i,a5.z);}`;
let blit_vert = `attribute vec2 a1;void main(){gl_Position=vec4(a1,0,1);}`;
let main_frag = `precision highp float;varying vec3 af,a6;varying float a8;uniform sampler2D a4[7];uniform vec3 ag;void main(){vec3 v=floor(fract(af*.005+.01)*64.);vec2 a=(v.xy+vec2(mod(v.z,8.),floor(v.z/8.))*64.+vec2(.5))/512.;vec4 i;int t=int(a8+.5);if(t==0)i=texture2D(a4[0],a);if(t==1)i=texture2D(a4[1],a);if(t==2)i=texture2D(a4[2],a);if(t==3)i=texture2D(a4[3],a);if(t==4)i=texture2D(a4[4],a);if(t==5)i=texture2D(a4[5],a);if(t==6)i=texture2D(a4[6],a);float r=dot(a6,normalize(vec3(5,1,5)));vec3 c=mix(mix(vec3(.9,.42,.44)+.1,vec3(1),.3),mix(vec3(1,.86,.39),vec3(1),.3),.5+.5*r);float n=.75+.25*max(0.,r),f=length(af.xz-ag.xz);if(f<20.&&af.y<ag.y)n*=1.-clamp(.5-.1*(f-5.),0.,.4);gl_FragColor=vec4(i.xyz*c*n,1);}`;
let main_vert = `attribute vec3 a1,ah;attribute float a9;uniform mat4 a2,a3;varying vec3 af,a6;varying float a8;void main(){af=a1,a6=normalize((a3*vec4(ah,0)).xyz),a8=a9,gl_Position=a2*a3*vec4(a1,1);}`;
let sky_frag = `precision highp float;varying vec3 af;void main(){vec3 m=normalize(af);float l=dot(normalize(vec3(5,1,5)),m),g=max(0.,500.*l-495.);vec3 x=mix(vec3(.13,.16,.31),vec3(.9,.42,.44),1.-2.*m.y),u=mix(mix(x,vec3(1,.86,.39),clamp(exp(l-1.)-.5*m.y,0.,1.)),vec3(1),g);gl_FragColor=vec4(u,1);}`;
let sky_vert = `attribute vec3 a1,ah;attribute float a9;uniform mat4 aa;varying vec3 af;void main(){af=a1,gl_Position=(aa*vec4(a1,1)).xyww;}`;

// ------------------------------------------------------------------------------------
// Object models
let skyboxGeo = csgSolidBake(csgSolidBox(0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0))[0];
let playerGeo = csgSolidBake(csgSolidBox(3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50))[0];
//export let worldSourceList:[number,string[]][]=[[0,["line","2","400","0","0","800","80","80","90","90","0"]],[0,["box","2","0","0","0","60","200","200","0","0","0","0"]],[0,["sub"]],[0,["box","2","284","0","0","208","200","200","0","0","0","0"]],[0,["sub"]],[0,["box","2","-284","0","0","208","200","200","0","0","0","0"]],[0,["sub"]],[0,["line","4","0","0","0","160","76","48","0","-90","0"]],[0,["add"]],[0,["box","4","0","0","-360","200","200","200","0","0","0","0"]],[0,["sub"]],[0,["line","5","0","0","0","160","40","40","0","-90","0"]],[0,["sub"]],[0,[""]]]
let [cannonGeo, _unused0] = csgSolidBake(csgSolidOpSubtract(csgSolidOpSubtract(csgSolidOpUnion(csgSolidOpSubtract(csgSolidOpSubtract(csgSolidOpSubtract(csgSolidLine(2, 400, 0, 0, 800, 80, 80, 90, 90, 0), csgSolidBox(2, 0, 0, 0, 60, 200, 200, 0, 0, 0, 0)), csgSolidBox(2, 284, 0, 0, 208, 200, 200, 0, 0, 0, 0)), csgSolidBox(2, -284, 0, 0, 208, 200, 200, 0, 0, 0, 0)), csgSolidLine(4, 0, 0, 0, 160, 76, 48, 0, -90, 0)), csgSolidBox(4, 0, 0, -360, 200, 200, 200, 0, 0, 0, 0)), csgSolidLine(5, 0, 0, 0, 160, 40, 40, 0, -90, 0)));
//export let worldSourceList:[number,string[]][]=[[0,[""]],[0,["#","Main","body","and","roof"]],[0,["box","6","0","0","0","130","100","130","0","0","0","0"]],[0,["box","6","0","110","0","100","25","100","0","0","0","0"]],[0,["sub"]],[0,[""]],[0,["#","Turrets"]],[0,["line","6","100","-119","100","400","50","50","0","0","0"]],[0,["add"]],[0,["line","6","100","-120","-100","400","50","50","0","0","0"]],[0,["add"]],[0,["line","6","-100","-120","-100","400","50","50","0","0","0"]],[0,["add"]],[0,["line","6","-100","-120","100","400","50","50","0","0","0"]],[0,["add"]],[0,[""]],[0,["#","Chop","the","bottom","off"]],[0,["box","6","0","-150","0","200","100","200","0","0","0","0"]],[0,["sub"]],[0,[""]],[0,["#","Turret","interior"]],[0,["line","4","100","-120","100","340","40","40","0","0","0"]],[0,["box","4","100","-50","100","50","170","50","0","0","0","0"]],[0,["sub"]],[0,["sub"]],[0,["line","4","100","0","-100","340","40","40","0","0","0"]],[0,["box","4","100","-50","-100","50","170","50","0","0","0","0"]],[0,["sub"]],[0,["sub"]],[0,["line","4","-100","-120","-100","340","40","40","0","0","0"]],[0,["box","4","-100","-50","-100","50","170","50","0","0","0","0"]],[0,["sub"]],[0,["sub"]],[0,["line","4","-100","-120","100","340","40","40","0","0","0"]],[0,["box","4","-100","-50","100","50","170","50","0","0","0","0"]],[0,["sub"]],[0,["sub"]],[0,[""]],[0,["#","Chop","the","top","off"]],[0,["box","6","0","250","0","200","100","200","0","0","0","0"]],[0,["sub"]],[0,[""]],[0,["#","Cut","doors","and","interior"]],[0,["box","5","0","-75","0","1","50","200","0","0","0","30"]],[0,["sub"]],[0,["box","5","0","-75","0","1","50","200","90","0","0","30"]],[0,["sub"]],[0,["box","5","0","-30","0","100","100","100","0","0","0","0"]],[0,["sub"]],[0,[""]],[0,["#","Turret","notches"]],[0,["box","4","0","150","100","200","20","10","0","0","0","0"]],[0,["sub"]],[0,["box","4","0","150","-100","200","20","10","0","0","0","0"]],[0,["sub"]],[0,["box","4","100","150","0","10","20","200","0","0","0","0"]],[0,["sub"]],[0,["box","4","-100","150","0","10","20","200","0","0","0","0"]],[0,["sub"]],[0,[""]],[0,[""]]]
let castleBase = csgSolidOpSubtract(csgSolidOpSubtract(csgSolidOpSubtract(csgSolidOpSubtract(csgSolidOpSubtract(csgSolidOpSubtract(csgSolidOpSubtract(csgSolidOpSubtract(csgSolidOpSubtract(csgSolidOpSubtract(csgSolidOpSubtract(csgSolidOpSubtract(csgSolidOpSubtract(csgSolidOpUnion(csgSolidOpUnion(csgSolidOpUnion(csgSolidOpUnion(csgSolidOpSubtract(csgSolidBox(6, 0, 0, 0, 130, 100, 130, 0, 0, 0, 0), csgSolidBox(6, 0, 110, 0, 100, 25, 100, 0, 0, 0, 0)), csgSolidLine(6, 100, -119, 100, 400, 50, 50, 0, 0, 0)), csgSolidLine(6, 100, -120, -100, 400, 50, 50, 0, 0, 0)), csgSolidLine(6, -100, -120, -100, 400, 50, 50, 0, 0, 0)), csgSolidLine(6, -100, -120, 100, 400, 50, 50, 0, 0, 0)), csgSolidBox(6, 0, -150, 0, 200, 100, 200, 0, 0, 0, 0)), csgSolidOpSubtract(csgSolidLine(4, 100, -120, 100, 340, 40, 40, 0, 0, 0), csgSolidBox(4, 100, -50, 100, 50, 170, 50, 0, 0, 0, 0))), csgSolidOpSubtract(csgSolidLine(4, 100, 0, -100, 340, 40, 40, 0, 0, 0), csgSolidBox(4, 100, -50, -100, 50, 170, 50, 0, 0, 0, 0))), csgSolidOpSubtract(csgSolidLine(4, -100, -120, -100, 340, 40, 40, 0, 0, 0), csgSolidBox(4, -100, -50, -100, 50, 170, 50, 0, 0, 0, 0))), csgSolidOpSubtract(csgSolidLine(4, -100, -120, 100, 340, 40, 40, 0, 0, 0), csgSolidBox(4, -100, -50, 100, 50, 170, 50, 0, 0, 0, 0))), csgSolidBox(6, 0, 250, 0, 200, 100, 200, 0, 0, 0, 0)), csgSolidBox(5, 0, -75, 0, 1, 50, 200, 0, 0, 0, 30)), csgSolidBox(5, 0, -75, 0, 1, 50, 200, 90, 0, 0, 30)), csgSolidBox(5, 0, -30, 0, 100, 100, 100, 0, 0, 0, 0)), csgSolidBox(4, 0, 150, 100, 200, 20, 10, 0, 0, 0, 0)), csgSolidBox(4, 0, 150, -100, 200, 20, 10, 0, 0, 0, 0)), csgSolidBox(4, 100, 150, 0, 10, 20, 200, 0, 0, 0, 0)), csgSolidBox(4, -100, 150, 0, 10, 20, 200, 0, 0, 0, 0));
let [castleGeo, _unused1] = csgSolidBake(castleBase);
let castleGibs = Array(8);
{
    const EXPLODE_CASTLE_TEXTURE = 1;
    let i = 0;
    for (let x = -100; x < 150; x += 200) {
        for (let y = -50; y < 200; y += 200) {
            for (let z = -100; z < 150; z += 200) {
                let build = csgSolidBake(csgSolidOpIntersect(castleBase, csgSolidBox(EXPLODE_CASTLE_TEXTURE, x, y, z, 80, 80, 80, Math.random() * 360, Math.random() * 360, Math.random() * 360, 20)));
                castleGibs[i++] = build[0];
            }
        }
    }
}
let gibCastle = (pos, vel) => {
    let ret = [];
    let i = 0;
    for (let x = -100; x < 150; x += 200) {
        for (let y = -50; y < 200; y += 200) {
            for (let z = -100; z < 150; z += 200) {
                ret.push({
                    kind: i++,
                    pos: v3AddScale(pos, [x, y, z], 0.25),
                    vel: v3AddScale(v3AddScale([0, 0, 0], vel, 0.8 / 33), [x, (y + 100) / 2, z], 1e-3),
                    offset: v3Negate([x, y, z]),
                    axis: v3Normalize([
                        2 * Math.random() - 1,
                        2 * Math.random() - 1,
                        2 * Math.random() - 1,
                    ]),
                    omega: 0.01 * Math.random(),
                    rotation: m4Ident,
                });
            }
        }
    }
    return ret;
};
// Hard qp level
// export let worldSourceList:[number,string[]][]=[[0,["box","0","0","-1000","1500","1026","1000","2000","0","0","0","0"]],[0,[""]],[0,["box","1","0","321","1211","1326","360","835","0","0","0","200"]],[0,["sub"]],[0,[""]],[0,[""]],[0,[""]],[0,["castle","-200","300","2200"]],[0,["castle","200","300","2200"]],[0,["castle","0","40","3000"]],[0,[""]]]
// let loadLevel4=()=>[csgSolidBake(csgSolidOpSubtract(csgSolidBox(0,0,-1000,1500,1026,1000,2000,0,0,0,0),csgSolidBox(1,0,321,1211,1326,360,835,0,0,0,200))),[[-200,300,2200],[200,300,2200],[0,40,3000]]]
// ------------------------------------------------------------------------------------
// World models
//export let worldSourceList:[number,string[]][]=[[0,["box","0","0","-1283","0","20000","1500","20000","0","0","0","0"]],[0,["box","0","0","1200","0","200","1000","1000","0","0","0","200"]],[0,["sub"]],[0,[""]],[0,[""]],[0,["castle","0","30","1000"]],[0,[""]],[0,[""]]]
let loadLevelBegin1 = () => [csgSolidBake(csgSolidOpSubtract(csgSolidBox(0, 0, -1283, 0, 20000, 1500, 20000, 0, 0, 0, 0), csgSolidBox(0, 0, 1200, 0, 200, 1000, 1000, 0, 0, 0, 200))), [[0, 30, 1000]]];
//export let worldSourceList:[number,string[]][]=[[0,["box","0","0","-891","1498","500","1000","2000","0","0","0","0"]],[0,["box","1","0","500","317","150","300","900","0","0","0","200"]],[0,["sub"]],[0,["box","2","0","-461","3875","450","1000","500","0","0","0","0"]],[0,["add"]],[0,[""]],[0,["castle","0","30","1000"]],[0,["#castle","0","120","2000"]],[0,["castle","0","140","3000"]],[0,[""]]]
let loadLevelBegin2 = () => [csgSolidBake(csgSolidOpUnion(csgSolidOpSubtract(csgSolidBox(0, 0, -891, 1498, 500, 1000, 2000, 0, 0, 0, 0), csgSolidBox(1, 0, 500, 317, 150, 300, 900, 0, 0, 0, 200)), csgSolidBox(2, 0, -461, 3875, 450, 1000, 500, 0, 0, 0, 0))), [[0, 30, 1000], [0, 140, 3000]]];
//export let worldSourceList:[number,string[]][]=[[0,["box","0","0","-903","1498","1026","1000","2000","0","0","0","0"]],[0,["box","1","0","500","317","150","300","900","0","0","0","200"]],[0,["sub"]],[0,["#box","2","0","-461","2375","1000","1000","1000","0","0","0","0"]],[0,["#box","2","0","-461","2175","900","1100","1100","0","0","0","0"]],[0,["#sub"]],[0,["#add"]],[0,["box","2","950","-461","2375","50","1000","1000","0","0","0","0"]],[0,["add"]],[0,["box","2","-950","-461","2375","50","1000","1000","0","0","0","0"]],[0,["add"]],[0,["box","2","0","-461","3342","50","1000","1000","90","0","0","0"]],[0,["add"]],[0,[""]],[0,["castle","500","140","3000"]],[0,["castle","-500","140","3000"]],[0,["castle","0","140","3000"]],[0,[""]]]
let loadLevelEasy3 = () => [csgSolidBake(csgSolidOpUnion(csgSolidOpUnion(csgSolidOpUnion(csgSolidOpSubtract(csgSolidBox(0, 0, -903, 1498, 1026, 1000, 2000, 0, 0, 0, 0), csgSolidBox(1, 0, 500, 317, 150, 300, 900, 0, 0, 0, 200)), csgSolidBox(2, 950, -461, 2375, 50, 1000, 1000, 0, 0, 0, 0)), csgSolidBox(2, -950, -461, 2375, 50, 1000, 1000, 0, 0, 0, 0)), csgSolidBox(2, 0, -461, 3342, 50, 1000, 1000, 90, 0, 0, 0))), [[500, 140, 3000], [-500, 140, 3000], [0, 140, 3000]]];
//export let worldSourceList:[number,string[]][]=[[0,["box","0","0","-1000","598","300","1000","915","0","0","0","0"]],[0,["box","0","-715","-1000","1213","300","1000","584","90","0","0","0"]],[0,["add"]],[0,["box","0","-1413","-1000","598","300","1000","915","0","0","0","0"]],[0,["add"]],[0,[""]],[0,["box","2","-1684","-238","1406","382","360","271","48","0","0","0"]],[0,["add"]],[0,[""]],[0,[""]],[0,["castle","0","38","1213"]],[0,["castle","-1415","38","1213"]],[0,["castle","-1415","38","0"]],[0,[""]]]
let loadLevelBraking = () => [csgSolidBake(csgSolidOpUnion(csgSolidOpUnion(csgSolidOpUnion(csgSolidBox(0, 0, -1000, 598, 300, 1000, 915, 0, 0, 0, 0), csgSolidBox(0, -715, -1000, 1213, 300, 1000, 584, 90, 0, 0, 0)), csgSolidBox(0, -1413, -1000, 598, 300, 1000, 915, 0, 0, 0, 0)), csgSolidBox(2, -1684, -238, 1406, 382, 360, 271, 48, 0, 0, 0))), [[0, 38, 1213], [-1415, 38, 1213], [-1415, 38, 0]]];
//export let worldSourceList:[number,string[]][]=[[0,["box","0","0","-1000","1500","1026","1000","2000","0","0","0","0"]],[0,[""]],[0,["box","1","0","85","978","328","360","494","0","0","0","0"]],[0,["add"]],[0,["box","1","0","85","1787","328","1154","494","0","0","0","0"]],[0,["add"]],[0,["box","1","0","799","-66","328","20","494","0","0","0","0"]],[0,["add"]],[0,["box","1","0","1185","-348","328","385","212","0","0","0","0"]],[0,["add"]],[0,[""]],[0,["castle","0","950","500"]],[0,["castle","0","600","500"]],[0,["castle","0","1250","500"]]]
let loadLevelStair = () => [csgSolidBake(csgSolidOpUnion(csgSolidOpUnion(csgSolidOpUnion(csgSolidOpUnion(csgSolidBox(0, 0, -1000, 1500, 1026, 1000, 2000, 0, 0, 0, 0), csgSolidBox(1, 0, 85, 978, 328, 360, 494, 0, 0, 0, 0)), csgSolidBox(1, 0, 85, 1787, 328, 1154, 494, 0, 0, 0, 0)), csgSolidBox(1, 0, 799, -66, 328, 20, 494, 0, 0, 0, 0)), csgSolidBox(1, 0, 1185, -348, 328, 385, 212, 0, 0, 0, 0))), [[0, 950, 500], [0, 600, 500], [0, 1250, 500]]];
//export let worldSourceList:[number,string[]][]=[[0,["box","0","0","-800","1100","1026","1000","1451","0","0","0","0"]],[0,[""]],[0,["box","1","0","521","-230","1326","360","1540","0","0","0","200"]],[0,["sub"]],[0,[""]],[0,[""]],[0,[""]],[0,[""]],[0,["castle","-400","450","1500"]],[0,["castle","0","500","1500"]],[0,["castle","400","450","1500"]],[0,[""]],[0,[""]]]
let loadLevelQuarter = () => [csgSolidBake(csgSolidOpSubtract(csgSolidBox(0, 0, -800, 1100, 1026, 1000, 1451, 0, 0, 0, 0), csgSolidBox(1, 0, 521, -230, 1326, 360, 1540, 0, 0, 0, 200))), [[-400, 450, 1500], [0, 500, 1500], [400, 450, 1500]]];
//export let worldSourceList:[number,string[]][]=[[0,["box","0","0","-903","336","507","1000","497","0","0","0","0"]],[0,["box","1","0","500","-295","150","300","900","0","0","0","200"]],[0,["sub"]],[0,[""]],[0,["box","2","0","-550","1266","330","100","652","0","0","0","0"]],[0,["add"]],[0,[""]],[0,["box","0","0","-903","2147","507","1000","497","0","0","0","0"]],[0,["add"]],[0,["box","1","0","500","2762","150","300","900","0","0","0","200"]],[0,["sub"]],[0,["box","5","0","-200","3504","1","1","900","0","0","0","130"]],[0,["sub"]],[0,[""]],[0,[""]],[0,[""]],[0,["castle","0","-309","1248"]],[0,["castle","0","-275","2584"]],[0,[""]]]
let loadLevelBounce = () => [csgSolidBake(csgSolidOpSubtract(csgSolidOpSubtract(csgSolidOpUnion(csgSolidOpUnion(csgSolidOpSubtract(csgSolidBox(0, 0, -903, 336, 507, 1000, 497, 0, 0, 0, 0), csgSolidBox(1, 0, 500, -295, 150, 300, 900, 0, 0, 0, 200)), csgSolidBox(2, 0, -550, 1266, 330, 100, 652, 0, 0, 0, 0)), csgSolidBox(0, 0, -903, 2147, 507, 1000, 497, 0, 0, 0, 0)), csgSolidBox(1, 0, 500, 2762, 150, 300, 900, 0, 0, 0, 200)), csgSolidBox(5, 0, -200, 3504, 1, 1, 900, 0, 0, 0, 130))), [[0, -309, 1248], [0, -275, 2584]]];
//export let worldSourceList:[number,string[]][]=[[0,["box","0","37","-903","576","535","1000","906","0","0","0","0"]],[0,["box","1","0","463","-411","150","300","900","0","0","0","200"]],[0,["sub"]],[0,[""]],[0,["box","2","-488","92","576","39","1130","1081","0","0","0","0"]],[0,["add"]],[0,[""]],[0,["line","2","-264","476","1036","500","200","200","90","90","0"]],[0,["sub"]],[0,[""]],[0,["castle","0","320","830"]],[0,["castle","0","462","1030"]],[0,["castle","-1000","462","1030"]],[0,[""]]]
let loadLevelRing = () => [csgSolidBake(csgSolidOpSubtract(csgSolidOpUnion(csgSolidOpSubtract(csgSolidBox(0, 37, -903, 576, 535, 1000, 906, 0, 0, 0, 0), csgSolidBox(1, 0, 463, -411, 150, 300, 900, 0, 0, 0, 200)), csgSolidBox(2, -488, 92, 576, 39, 1130, 1081, 0, 0, 0, 0)), csgSolidLine(2, -264, 476, 1036, 500, 200, 200, 90, 90, 0))), [[0, 320, 830], [0, 462, 1030], [-1000, 462, 1030]]];
let loadLevelX = () => [csgSolidBake(csgSolidOpUnion(csgSolidOpSubtract(csgSolidBox(1, 0, 201, 1371, 1376, 2029, 1008, 0, 0, 0, 0), csgSolidBox(1, 0, 0, 579, 1, 1, 1, 0, 0, 0, 400)), csgSolidBox(2, 0, -50, 0, 50, 50, 50, 0, 0, 0, 0))), [[230, 201, 446], [-230, 201, 446], [210, -253, 446], [-210, -253, 446]]];
const START_LEVEL = 0;
let levelLoaders = [
    loadLevelBegin1,
    loadLevelBegin2,
    loadLevelEasy3,
    loadLevelBraking,
    loadLevelStair,
    loadLevelBounce,
    loadLevelRing,
    loadLevelQuarter,
    loadLevelX,
];
let lastLevel = levelLoaders.length - 1;
let worldCastles;
let worldGeo;
let worldFn;
let worldGetGeo = () => worldGeo;
let worldGetCastles = () => worldCastles;
let loadLevel = (i) => {
    [[worldGeo, worldFn], worldCastles] = levelLoaders[i]();
};
// ------------------------------------------------------------------------------------
let startingAmmos = [
    2,
    2,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
];
let hints = [
    'SMASH THE CASTLE',
    'SMASH ALL THE CASTLES',
    'RE-USE THE CANNON',
    'USE THE CANNON TO SLOW DOWN',
    'SHOOT INTO THE AIR',
    'MIND THE GAP',
    'BLAST THROUGH THE WALL',
    'USE THE RAMP',
    'RIDE THE CURVE',
];
const eps = 0.001;
let worldSampleNormal = (pos) => v3Normalize([
    worldFn(v3Add(pos, [eps, 0, 0]))[1] - worldFn(v3Sub(pos, [eps, 0, 0]))[1],
    worldFn(v3Add(pos, [0, eps, 0]))[1] - worldFn(v3Sub(pos, [0, eps, 0]))[1],
    worldFn(v3Add(pos, [0, 0, eps]))[1] - worldFn(v3Sub(pos, [0, 0, eps]))[1]
]);
let worldNearestSurfacePoint = (pos) => {
    for (let i = 0, marchPoint = pos, dist, tag, norm; i < 50; ++i) {
        [tag, dist] = worldFn(marchPoint);
        norm = worldSampleNormal(marchPoint);
        if (dist < eps) {
            dist = v3Sub(pos, marchPoint);
            return [marchPoint, v3Dot2(dist) < eps ? norm : v3Normalize(dist), v3Length(dist), tag];
        }
        marchPoint = v3AddScale(marchPoint, norm, -dist);
    }
    return undefined;
};
let worldRaycast = (pos, normalizedDir, len) => {
    let i = 0, traveled = 0, marchPoint = pos, dist;
    for (; i < 50 && traveled < len; ++i) {
        dist = worldFn(marchPoint)[1];
        traveled += dist;
        if (dist < eps) {
            return traveled;
        }
        marchPoint = v3AddScale(marchPoint, normalizedDir, dist);
    }
    return len;
};

// --------------------------------------------------------
// ZzFXMicro - Zuper Zmall Zound Zynth - v1.1.8
// https://github.com/KilledByAPixel/ZzFX
let zzfxV = .3; // volume
let zzfxR = 44100; // sample rate
let zzfxX = new (window.AudioContext || webkitAudioContext); // audio context
let zzfxP = (samples) => // play samples
 {
    // create buffer and source
    let buffer = zzfxX.createBuffer(1, samples.length, zzfxR), source = zzfxX.createBufferSource();
    // copy samples to buffer and play
    buffer.getChannelData(0).set(samples);
    source.buffer = buffer;
    source.connect(zzfxX.destination);
    source.start();
    return source;
};
let zzfxG = // generate samples
 (
// parameters
volume = 1, randomness = .05, frequency = 220, attack = 0, sustain = 0, release = .1, shape = 0, shapeCurve = 1, slide = 0, deltaSlide = 0, pitchJump = 0, pitchJumpTime = 0, repeatTime = 0, noise = 0, modulation = 0, bitCrush = 0, delay = 0, sustainVolume = 1, decay = 0, tremolo = 0) => {
    // init parameters
    let PI2 = Math.PI * 2, sign = (v) => v > 0 ? 1 : -1, startSlide = slide *= 500 * PI2 / zzfxR / zzfxR, startFrequency = frequency *= (1 + randomness * 2 * Math.random() - randomness)
        * PI2 / zzfxR, b = [], t = 0, tm = 0, i = 0, j = 1, r = 0, c = 0, s = 0, f, length;
    // scale by sample rate
    attack = attack * zzfxR + 9; // minimum attack to prevent pop
    decay *= zzfxR;
    sustain *= zzfxR;
    release *= zzfxR;
    delay *= zzfxR;
    deltaSlide *= 500 * PI2 / zzfxR ** 3;
    modulation *= PI2 / zzfxR;
    pitchJump *= PI2 / zzfxR;
    pitchJumpTime *= zzfxR;
    repeatTime = repeatTime * zzfxR | 0;
    // generate waveform
    for (length = attack + decay + sustain + release + delay | 0; i < length; b[i++] = s) {
        if (!(++c % (bitCrush * 100 | 0))) // bit crush
         {
            s = shape ? shape > 1 ? shape > 2 ? shape > 3 ? // wave shape
                Math.sin((t % PI2) ** 3) : // 4 noise
                Math.max(Math.min(Math.tan(t), 1), -1) : // 3 tan
                1 - (2 * t / PI2 % 2 + 2) % 2 : // 2 saw
                1 - 4 * Math.abs(Math.round(t / PI2) - t / PI2) : // 1 triangle
                Math.sin(t); // 0 sin
            s = (repeatTime ?
                1 - tremolo + tremolo * Math.sin(PI2 * i / repeatTime) // tremolo
                : 1) *
                sign(s) * (Math.abs(s) ** shapeCurve) * // curve 0=square, 2=pointy
                volume * zzfxV * ( // envelope
            i < attack ? i / attack : // attack
                i < attack + decay ? // decay
                    1 - ((i - attack) / decay) * (1 - sustainVolume) : // decay falloff
                    i < attack + decay + sustain ? // sustain
                        sustainVolume : // sustain volume
                        i < length - delay ? // release
                            (length - i - delay) / release * // release falloff
                                sustainVolume : // release volume
                            0); // post release
            s = delay ? s / 2 + (delay > i ? 0 : // delay
                (i < length - delay ? 1 : (length - i) / delay) * // release delay
                    b[i - delay | 0] / 2) : s; // sample delay
        }
        f = (frequency += slide += deltaSlide) * // frequency
            Math.cos(modulation * tm++); // modulation
        t += f - f * noise * (1 - (Math.sin(i) + 1) * 1e9 % 2); // noise
        if (j && ++j > pitchJumpTime) // pitch jump
         {
            frequency += pitchJump; // apply pitch jump
            startFrequency += pitchJump; // also apply to start
            j = 0; // reset pitch jump time
        }
        if (repeatTime && !(++r % repeatTime)) // repeat
         {
            frequency = startFrequency; // reset frequency
            slide = startSlide; // reset slide
            j = j || 1; // reset pitch jump time
        }
    }
    return b;
};
// --------------------------------------------------------
// sound defs:
let sndGround = zzfxG(...[1.27, , 478, .02, .02, .06, 1, 2.75, , , , , , .9, 45, .5, , .59, , .26]);
let sndShoot = zzfxG(...[1.37, , 827, .02, .01, .32, 4, 4.22, , .5, , , , 1.7, , .6, , .47]);
let sndCastle = zzfxG(...[2, .1, 185, .03, , .61, 4, 4, , , , , , 1, , 1, .5, .8, .1, .1]);
let sndWin = zzfxG(...[2.11, 0, 97.99886, .08, .7, .34, , 1.92, -0.5, 10, 50, .02, .12, , , .1, .06, .17, .01, .66]);
let sndLose = zzfxG(...[2.11, 0, 261.6256, .08, .76, .35, , 1.92, , , -50, -0.01, .2, , , .1, .06, .27, .09, .2]);
let sndBird = sndWin; //zzfxG(...[,,146.8324,,.25,.38,,1.47,-3.9,,70,.03,.03,,5,.1,,.97,.23,.41])

let gameStateNew = () => ({
    mode: 0 /* Menu */,
    holdingMouse: False,
    lockView: False,
    yaw: Math.PI,
    pitch: 0,
    pos: [0, 10, 0],
    vel: [0, 0, 0],
    rotSpeed: 0,
    rotAxis: [0, 0, 0],
    camBack: 100,
    ammo: startingAmmos[START_LEVEL],
    castlesHit: [],
    modeTick: 0,
    level: START_LEVEL,
    ungrounded: 0,
});
let gameStateLerp = (a, b, t) => ({
    mode: b.mode,
    holdingMouse: b.holdingMouse,
    lockView: b.lockView,
    yaw: b.yaw,
    pitch: b.pitch,
    pos: vecLerp(a.pos, b.pos, t),
    vel: b.vel,
    rotSpeed: b.rotSpeed,
    rotAxis: b.rotAxis,
    camBack: lerp(a.camBack, b.camBack, t),
    ammo: b.ammo,
    castlesHit: b.castlesHit,
    modeTick: b.modeTick,
    level: b.level,
    ungrounded: b.ungrounded,
});
let predictShot = (yaw, pitch, pos, castlesHit) => {
    let lookVec = m4MulPoint(m4Mul(m4RotY(yaw), m4RotX(-pitch)), [0, 0, -1]);
    let ret = new Float32Array(6 * 150);
    let loopout;
    let castles = worldGetCastles();
    let vel = v3AddScale([0, 0, 0], lookVec, 30);
    for (let i = 0, j = 0; i < 150; ++i) {
        ret[j++] = pos[0];
        ret[j++] = pos[1];
        ret[j++] = pos[2];
        if (!loopout) {
            vel = v3Add(vel, [0, -0.6, 0]);
            pos = v3Add(pos, vel);
            let testpos = pos;
            for (let i = 0; i < 4; ++i) {
                testpos = v3AddScale(testpos, vel, 0.25);
                for (let i = 0; i < castles.length; ++i) {
                    if (castlesHit.indexOf(i) >= 0)
                        continue;
                    if (v3Dot2(v3Sub(castles[i], testpos)) < 3600) {
                        loopout = 1;
                    }
                }
            }
            let [nearPos, nearNorm, nearDist] = worldNearestSurfacePoint(pos);
            if (nearDist < 10 && v3Dot(nearNorm, vel) < 0) {
                pos = v3AddScale(nearPos, nearNorm, 10);
                loopout = 1;
            }
        }
        ret[j++] = pos[0];
        ret[j++] = pos[1];
        ret[j++] = pos[2];
    }
    return [ret, pos];
};
let firstLaterAim = true;
let gameStateTick = (prevState, inputs) => {
    let state = gameStateLerp(prevState, prevState, 0);
    state.modeTick++;
    state.ungrounded++;
    if (state.mode != 0 /* Menu */ && state.mode != 4 /* Dead */ && state.mode != 5 /* Win */) {
        // state.lockView = (inputs.keysDown[2] && (state.mode == GameMode.FirstAim || state.mode == GameMode.LaterAim)) as any
        state.yaw += inputs.mouseAccX * 0.005;
        state.pitch += inputs.mouseAccY * 0.005;
        state.pitch = Math.max(-1.57, Math.min(1.57, state.pitch));
        if (state.yaw < 0)
            state.yaw += 2 * Math.PI;
        if (state.yaw > 2 * Math.PI)
            state.yaw -= 2 * Math.PI;
    }
    let lookVec = m4MulPoint(m4Mul(m4RotY(state.yaw), m4RotX(-state.pitch)), [0, 0, -1]);
    let click = !state.holdingMouse && inputs.keysDown[0];
    state.holdingMouse = inputs.keysDown[0];
    let castles = worldGetCastles();
    if (state.mode == 3 /* Ball */ && (inputs.keysDown['R'] || state.pos[1] < -1000)) {
        state.mode = 4 /* Dead */;
        state.modeTick = 0;
        state.ammo = 0;
        return state;
    }
    if (state.mode == 0 /* Menu */) {
        if (clickedIn.a) {
            zzfxP(sndWin);
            state.mode = 2 /* LaterAim */;
        }
    }
    else if (state.mode == 1 /* FirstAim */ || state.mode == 2 /* LaterAim */) {
        state.vel = v3AddScale(state.vel, [0, -0.6, 0], 0.5);
        let testpos = state.pos;
        for (let i = 0; i < 4; ++i) {
            testpos = v3AddScale(testpos, state.vel, 0.25);
            for (let i = 0; i < castles.length; ++i) {
                if (state.castlesHit.indexOf(i) >= 0)
                    continue;
                if (v3Dot2(v3Sub(castles[i], testpos)) < 3600) {
                    zzfxP(sndCastle);
                    state.castlesHit.push(i);
                }
            }
        }
        if (state.castlesHit.length == castles.length) {
            state.modeTick = 0;
            state.mode = 5 /* Win */;
        }
        if (((state.mode == 1 /* FirstAim */ || firstLaterAim) && click || (state.mode == 2 /* LaterAim */ && !firstLaterAim && !inputs.keysDown[0]) && state.ammo > 0)) {
            zzfxP(sndShoot);
            state.ammo -= 1;
            state.mode = 3 /* Ball */;
            state.vel = v3AddScale([0, 0, 0], lookVec, 30);
            state.rotSpeed = 0;
            firstLaterAim = false;
        }
        state.pos = v3AddScale(state.pos, state.vel, 0.5);
        let [nearPos, nearNorm, nearDist] = worldNearestSurfacePoint(state.pos);
        if (nearDist < 10 && v3Dot(nearNorm, state.vel) < 0) {
            state.pos = v3AddScale(nearPos, nearNorm, 10);
            state.vel = v3Reflect(state.vel, nearNorm, 0.2, 0.8);
        }
    }
    else if (state.mode == 3 /* Ball */ || state.mode == 4 /* Dead */ || state.mode == 5 /* Win */) {
        state.vel = v3Add(state.vel, [0, -0.6, 0]);
        if (state.mode == 3 /* Ball */) {
            let testpos = state.pos;
            for (let i = 0; i < 4; ++i) {
                testpos = v3AddScale(testpos, state.vel, 0.25);
                for (let i = 0; i < castles.length; ++i) {
                    if (state.castlesHit.indexOf(i) >= 0)
                        continue;
                    if (v3Dot2(v3Sub(castles[i], testpos)) < 3600) {
                        zzfxP(sndCastle);
                        state.castlesHit.push(i);
                    }
                }
            }
            if (state.castlesHit.length == castles.length) {
                state.modeTick = 0;
                state.mode = 5 /* Win */;
            }
            if (click) {
                state.modeTick = 0;
                state.mode = state.ammo > 0
                    ? 2 /* LaterAim */
                    : 4 /* Dead */;
            }
        }
        else {
            if (state.modeTick > 30) {
                if (state.mode == 5 /* Win */) {
                    if (state.level < lastLevel) {
                        if (state.ammo > 0) {
                            zzfxP(sndBird);
                        }
                        else {
                            zzfxP(sndWin);
                        }
                        let level = state.level + 1;
                        state = gameStateNew();
                        state.mode = 1 /* FirstAim */;
                        state.level = level;
                        state.ammo = startingAmmos[level];
                    }
                }
                else {
                    zzfxP(sndLose);
                    let level = state.level;
                    state = gameStateNew();
                    state.mode = 1 /* FirstAim */;
                    state.level = level;
                    state.ammo = startingAmmos[level];
                }
            }
        }
        state.pos = v3Add(state.pos, state.vel);
        let [nearPos, nearNorm, nearDist, tag] = worldNearestSurfacePoint(state.pos);
        if (nearDist < 10 && v3Dot(nearNorm, state.vel) < 0) {
            if (state.ungrounded > 10) {
                zzfxP(sndGround);
            }
            state.ungrounded = 0;
            state.pos = v3AddScale(nearPos, nearNorm, 10);
            state.vel = v3Reflect(state.vel, nearNorm, tag == 2 ? 0.9 : 0.2, 0.998);
            state.rotSpeed = v3Length(state.vel) / 10 / 33;
            state.rotAxis = v3Negate(v3Normalize(v3Cross(state.vel, nearNorm)));
        }
    }
    state.camBack = lerp(state.camBack, worldRaycast(state.pos, v3Negate(lookVec), 100), 0.5);
    return state;
};

let textures = Array(7).fill(0);
let genGrass = () => (_x, _y, _z, i, out) => {
    let brightness = Math.pow(Math.random(), 4.0);
    let darkness = Math.pow(Math.random(), 4.0);
    let adjust = 0.6 - 0.1 * brightness + 0.1 * darkness;
    out[i + 0] = 255 * 0.3;
    out[i + 1] = 255 * adjust;
    out[i + 2] = 255 * 0.25 * adjust;
};
let genBricks = () => (x, y, z, i, out) => {
    x += (2 * Math.random() - 1) / 64;
    y += (2 * Math.random() - 1) / 64;
    z += (2 * Math.random() - 1) / 64;
    let mortar = 0.05;
    x += 0.03;
    z += 0.03;
    x *= 4;
    z *= 4;
    y *= 8;
    x %= 1;
    z %= 1;
    y %= 1;
    let inner = x > mortar && x < (1 - mortar) &&
        z > mortar && z < (1 - mortar) &&
        y > (2 * mortar) && y < (1 - 2 * mortar);
    let brightness = 0.9;
    brightness -= brightness * (0.1 * Math.pow(Math.random(), 8));
    brightness += brightness * (0.1 * Math.pow(Math.random(), 8));
    if (inner) {
        out[i + 0] = 80 * brightness;
        out[i + 1] = 90 * brightness;
        out[i + 2] = 100 * brightness;
    }
    else {
        out[i + 0] = 80 * brightness;
        out[i + 1] = 70 * brightness;
        out[i + 2] = 60 * brightness;
    }
};
let genWood = () => {
    let noise = (v) => 0.25 * (Math.sin(10. * v)
        + .5 * Math.sin(17. * v)
        + .25 * Math.sin(21. * v)
        + .25 * Math.sin(23. * v)
        + Math.sin(29. * v)
        + .5 * Math.sin(31. * v)
        + .25 * Math.sin(37. * v)
        + .25 * Math.sin(51. * v));
    return (x, y, z, i, out) => {
        let darkness = Math.pow(Math.random(), 8.0);
        let b = 0.3 * (noise(5 * x) + noise(5 * z) + noise(0.2 * y)) * 0.25 + 0.75;
        b -= 0.05 * darkness;
        out[i + 0] = 255 * b;
        out[i + 1] = 220 * b;
        out[i + 2] = 200 * b;
    };
};
let genCobble = () => {
    let tilingDistance3D = (x1, y1, z1, x2, y2, z2) => (x1 = Math.abs(x1 - x2),
        x1 = Math.min(x1, 1 - x1),
        y1 = Math.abs(y1 - y2),
        y1 = Math.min(y1, 1 - y1),
        z1 = Math.abs(z1 - z2),
        z1 = Math.min(z1, 1 - z1),
        Math.hypot(x1, y1, z1));
    let voronoiCenters = Array(40).fill(undefined).map(_ => [Math.random(), Math.random(), Math.random()]);
    return (x, y, z, i, out) => {
        let minDistance = Infinity;
        let secondMinDistance = Infinity;
        for (let [cx, cy, cz] of voronoiCenters) {
            let dist = tilingDistance3D(x, y, z, cx, cy, cz);
            if (dist < minDistance) {
                secondMinDistance = minDistance;
                minDistance = dist;
            }
            else if (dist < secondMinDistance) {
                secondMinDistance = dist;
            }
        }
        let edgeDistance = Math.abs(minDistance - secondMinDistance);
        let brightness = edgeDistance < 0.04 ? 0.3 + 3 * edgeDistance : 0.4 + 0.6 * edgeDistance;
        brightness -= brightness * (0.1 * Math.pow(Math.random(), 8));
        brightness += brightness * (0.1 * Math.pow(Math.random(), 8));
        out[i + 0] = 230 * brightness;
        out[i + 1] = 200 * brightness;
        out[i + 2] = 255 * brightness;
    };
};
let genBallTex = () => (x, y, z, i, out) => {
    let brightness = Math.pow(Math.random(), 4.0);
    let darkness = Math.pow(Math.random(), 4.0);
    if (x > 0.5)
        x = 1.0 - x;
    if (y > 0.5)
        y = 1.0 - y;
    if (z > 0.5)
        z = 1.0 - z;
    let a = Math.hypot(x, y) < 0.05;
    let b = Math.hypot(x, z) < 0.05;
    let c = Math.hypot(y, z) < 0.05;
    let t = a || b || c ? 0.5 : 1;
    out[i + 0] = t * 100 * (0.6 - 0.1 * brightness + 0.1 * darkness);
    out[i + 1] = t * 150 * (0.6 - 0.1 * brightness + 0.1 * darkness);
    out[i + 2] = t * 180 * (0.6 - 0.1 * brightness + 0.1 * darkness);
};
let genCannonTex = () => (_x, _y, _z, i, out) => {
    let brightness = Math.pow(Math.random(), 4.0);
    let darkness = Math.pow(Math.random(), 4.0);
    let t = 1;
    out[i + 0] = t * 100 * (0.6 - 0.1 * brightness + 0.1 * darkness);
    out[i + 1] = t * 150 * (0.6 - 0.1 * brightness + 0.1 * darkness);
    out[i + 2] = t * 180 * (0.6 - 0.1 * brightness + 0.1 * darkness);
};
let genCannonTexDark = () => (_x, _y, _z, i, out) => {
    let brightness = Math.pow(Math.random(), 4.0);
    let darkness = Math.pow(Math.random(), 4.0);
    let t = 1;
    out[i + 0] = 0.5 * t * 100 * (0.6 - 0.1 * brightness + 0.1 * darkness);
    out[i + 1] = 0.5 * t * 150 * (0.6 - 0.1 * brightness + 0.1 * darkness);
    out[i + 2] = 0.5 * t * 180 * (0.6 - 0.1 * brightness + 0.1 * darkness);
};
let doGen = (fnBuilder) => {
    let fn = fnBuilder();
    let ret = new Uint8Array(512 * 512 * 4);
    for (let x = 0; x < 64; ++x) {
        for (let y = 0; y < 64; ++y) {
            for (let z = 0; z < 64; ++z) {
                let tx = z % 8;
                let ty = Math.floor(z / 8);
                let idx = 4 * (x + 64 * tx + 512 * (y + 64 * ty));
                fn(x / 64, y / 64, z / 64, idx, ret);
            }
        }
    }
    return ret;
};
let genTex3d = (texIdx, getColor) => {
    let code = `let a=(${doGen.toString()})(${getColor.toString()});self.postMessage(a,[a.buffer])`;
    let worker = new Worker(URL.createObjectURL(new Blob([code])));
    new Promise(resolve => {
        worker.onmessage = e => resolve(e.data);
    }).then((ret) => {
        let tex = G.createTexture();
        G.bindTexture(TEXTURE_2D, tex);
        G.texImage2D(TEXTURE_2D, 0, RGBA, 512, 512, 0, RGBA, UNSIGNED_BYTE, ret);
        G.texParameteri(TEXTURE_2D, TEXTURE_MIN_FILTER, NEAREST);
        G.texParameteri(TEXTURE_2D, TEXTURE_MAG_FILTER, NEAREST);
        G.texParameteri(TEXTURE_2D, TEXTURE_WRAP_S, REPEAT);
        G.texParameteri(TEXTURE_2D, TEXTURE_WRAP_T, REPEAT);
        textures[texIdx] = tex;
    });
};
let ret = new Uint8Array(512 * 512 * 4);
for (let i = 0; i < ret.length; ++i) {
    ret[i] = 255;
}
let fallback = G.createTexture();
G.bindTexture(TEXTURE_2D, fallback);
G.texImage2D(TEXTURE_2D, 0, RGBA, 512, 512, 0, RGBA, UNSIGNED_BYTE, ret);
G.texParameteri(TEXTURE_2D, TEXTURE_MIN_FILTER, NEAREST);
G.texParameteri(TEXTURE_2D, TEXTURE_MAG_FILTER, NEAREST);
G.texParameteri(TEXTURE_2D, TEXTURE_WRAP_S, REPEAT);
G.texParameteri(TEXTURE_2D, TEXTURE_WRAP_T, REPEAT);
genTex3d(0, genGrass);
genTex3d(1, genCobble);
genTex3d(2, genWood);
genTex3d(3, genBallTex);
genTex3d(4, genCannonTex);
genTex3d(5, genCannonTexDark);
genTex3d(6, genBricks);
let bindTextureUniforms = (shader) => {
    G.uniform1iv(G.getUniformLocation(shader, 'a4'), textures.map((tex, i) => (G.activeTexture(TEXTURE0 + i),
        G.bindTexture(TEXTURE_2D, tex || fallback),
        i)));
};

let ctx = C2.getContext('2d');
let ctxtx = G.createTexture();
let fade = 0;
G.bindTexture(TEXTURE_2D, ctxtx);
G.texParameteri(TEXTURE_2D, TEXTURE_MIN_FILTER, NEAREST);
G.texParameteri(TEXTURE_2D, TEXTURE_MAG_FILTER, NEAREST);
G.texParameteri(TEXTURE_2D, TEXTURE_WRAP_S, CLAMP_TO_EDGE);
G.texParameteri(TEXTURE_2D, TEXTURE_WRAP_T, CLAMP_TO_EDGE);
let blitTriBuffer = G.createBuffer();
G.bindBuffer(ARRAY_BUFFER, blitTriBuffer);
G.bufferData(ARRAY_BUFFER, Float32Array.of(-1, -1, 1, -1, -1, 1, 1, 1, -1, 1, 1, -1), STATIC_DRAW);
G.enable(DEPTH_TEST);
G.enable(BLEND);
G.blendFunc(SRC_ALPHA, ONE_MINUS_SRC_ALPHA);
G.depthFunc(LEQUAL);
let drawRayIndex = G.createBuffer();
G.bindBuffer(ELEMENT_ARRAY_BUFFER, drawRayIndex);
G.bufferData(ELEMENT_ARRAY_BUFFER, new Uint16Array(Array(2 * 150).fill(0).map((_, i) => i)), STATIC_DRAW);
let drawRayVertex = G.createBuffer();
G.bindBuffer(ARRAY_BUFFER, drawRayVertex);
G.bufferData(ARRAY_BUFFER, new Float32Array(6 * 150), DYNAMIC_DRAW);
let shaderCompile = (vert, frag) => {
    let vs = G.createShader(VERTEX_SHADER);
    let fs = G.createShader(FRAGMENT_SHADER);
    let shader = G.createProgram();
    G.shaderSource(vs, vert);
    G.compileShader(vs);
    G.shaderSource(fs, frag);
    G.compileShader(fs);
    G.attachShader(shader, vs);
    G.attachShader(shader, fs);
    G.linkProgram(shader);
    return shader;
};
let mainShader = shaderCompile(main_vert, main_frag);
let skyShader = shaderCompile(sky_vert, sky_frag);
let aimRayShader = shaderCompile(aimRay_vert, aimRay_frag);
let blitShader = shaderCompile(blit_vert, blit_frag);
let modelGeoDraw = (me, shaderProg) => {
    G.bindBuffer(ARRAY_BUFFER, me.vertexBuffer);
    let posLoc = G.getAttribLocation(shaderProg, 'a1');
    G.enableVertexAttribArray(posLoc);
    G.vertexAttribPointer(posLoc, 3, FLOAT, false, 0, 0);
    G.bindBuffer(ARRAY_BUFFER, me.normalBuffer);
    posLoc = G.getAttribLocation(shaderProg, 'ah');
    if (posLoc >= 0) {
        G.enableVertexAttribArray(posLoc);
        G.vertexAttribPointer(posLoc, 3, FLOAT, false, 0, 0);
        G.bindBuffer(ARRAY_BUFFER, me.tagBuffer);
        posLoc = G.getAttribLocation(shaderProg, 'a9');
        G.enableVertexAttribArray(posLoc);
        G.vertexAttribPointer(posLoc, 1, FLOAT, false, 0, 0);
    }
    G.bindBuffer(ELEMENT_ARRAY_BUFFER, me.indexBuffer);
    G.drawElements(TRIANGLES, me.indexBufferLen, UNSIGNED_SHORT, 0);
};
let ballRot = m4Ident;
let mainPitch = 0;
let mainYaw = 0;
let unlockT = 0;
let lazyPitch = 0;
let lazyYaw = 0;
let camOff = [0, 0, 0];
let cannonPos = [0, 0, 0];
let ballPos = [0, 0, 0];
let oldMode;
let modeT = 0;
let uihash;
let gibbedCastles = [];
let castleGibStates = [];
let resize = () => {
    G.viewport(0, 0, C2.width = CC.width = window.innerWidth / 2, C2.height = CC.height = window.innerHeight / 2);
    uihash = 9;
};
let renderGame = (earlyInputs, state, dt) => {
    if (state.mode == 4 /* Dead */ || state.mode == 5 /* Win */) {
        earlyInputs.mouseAccX = earlyInputs.mouseAccY = 0;
    }
    let predictedYaw = earlyInputs.mouseAccX * 0.005 + state.yaw;
    let predictedPitch = earlyInputs.mouseAccY * 0.005 + state.pitch;
    predictedPitch = Math.max(-1.5, Math.min(1.5, predictedPitch));
    modeT += dt;
    if (state.lockView) {
        unlockT = 500;
    }
    else if (unlockT > 0) {
        mainPitch = lerp(mainPitch, predictedPitch, 0.02 * dt);
        mainYaw = radLerp(mainYaw, predictedYaw, 0.02 * dt);
        unlockT -= dt;
    }
    else {
        mainPitch = predictedPitch;
        mainYaw = predictedYaw;
    }
    let lookVec = m4MulPoint(m4Mul(m4RotY(mainYaw), m4RotX(-mainPitch)), [0, 0, -state.camBack]);
    camOff = vecLerp(camOff, state.mode == 3 /* Ball */ || state.mode == 4 /* Dead */ || state.mode == 5 /* Win */ ? [0, 20, 0] : [-20, 30, 0], 0.01 * dt);
    ballRot = m4Mul(m4AxisAngle(state.rotAxis, state.rotSpeed * dt), ballRot);
    let lookMat = m4Mul(m4RotX(mainPitch), m4RotY(-mainYaw));
    let fwdLookMat = m4Mul(m4RotY(mainYaw), m4RotX(-mainPitch));
    let camOffset = m4MulPoint(fwdLookMat, camOff);
    let viewMat = m4Mul(lookMat, m4Translate(v3Sub(lookVec, v3Add(state.pos, camOffset))));
    let projectionMat = m4Perspective(CC.height / CC.width, 0.1, 10000);
    let vp = m4Mul(projectionMat, viewMat);
    let drawCannon = 0, modelMat;
    if (state.mode == 3 /* Ball */ || state.mode == 4 /* Dead */ || state.mode == 5 /* Win */) {
        ballPos = state.pos;
        drawCannon = (state.mode == 3 /* Ball */ && modeT < 100);
        drawBall(vp, mainShader, 0);
    }
    if (state.mode == 1 /* FirstAim */ || state.mode == 2 /* LaterAim */) {
        lazyPitch = lerp(lazyPitch, predictedPitch, 0.01 * dt);
        lazyYaw = radLerp(lazyYaw, predictedYaw, 0.01 * dt);
        cannonPos = state.pos;
        drawCannon = 1;
    }
    if (drawCannon || state.mode == 0 /* Menu */) {
        if (state.mode == 0 /* Menu */) {
            lazyYaw = 1;
            lazyPitch = -.5;
            cannonPos = state.pos;
        }
        modelMat = m4Mul(m4Mul(m4Translate(cannonPos), m4Mul(m4RotY(lazyYaw), m4RotX(-lazyPitch))), m4Scale(0.12));
        G.useProgram(mainShader);
        G.uniformMatrix4fv(G.getUniformLocation(mainShader, 'a3'), false, modelMat);
        G.uniformMatrix4fv(G.getUniformLocation(mainShader, 'a2'), false, vp);
        bindTextureUniforms(mainShader);
        modelGeoDraw(cannonGeo, mainShader);
    }
    let castles = worldGetCastles();
    for (let i = 0; i < castles.length; ++i) {
        if (state.castlesHit.indexOf(i) >= 0) {
            if (gibbedCastles.indexOf(i) < 0) {
                gibbedCastles.push(i);
                castleGibStates.push(...gibCastle(castles[i], state.vel));
            }
            continue;
        }
        modelMat = m4Mul(m4Mul(m4Translate(v3Add(castles[i], [0, 10 * Math.sin(modeT / 800), 0])), m4Scale(0.25)), m4Mul(m4RotY(modeT / 1000 + i), m4RotX(0.15)));
        G.useProgram(mainShader);
        G.uniformMatrix4fv(G.getUniformLocation(mainShader, 'a3'), false, modelMat);
        G.uniformMatrix4fv(G.getUniformLocation(mainShader, 'a2'), false, vp);
        bindTextureUniforms(mainShader);
        modelGeoDraw(castleGeo, mainShader);
    }
    for (let gib of castleGibStates) {
        let grav = -0.6 / 33 / 33;
        gib.rotation = m4Mul(m4AxisAngle(gib.axis, gib.omega * dt), gib.rotation);
        gib.vel = v3AddScale(gib.vel, [0, grav, 0], dt);
        gib.pos = v3AddScale(gib.pos, gib.vel, dt);
        modelMat = m4Mul(m4Mul(m4Translate(gib.pos), m4Scale(0.25)), m4Mul(gib.rotation, m4Translate(gib.offset)));
        G.useProgram(mainShader);
        G.uniformMatrix4fv(G.getUniformLocation(mainShader, 'a3'), false, modelMat);
        G.uniformMatrix4fv(G.getUniformLocation(mainShader, 'a2'), false, vp);
        bindTextureUniforms(mainShader);
        modelGeoDraw(castleGibs[gib.kind], mainShader);
    }
    if (oldMode != state.mode) {
        oldMode = state.mode;
        if (oldMode == 1 /* FirstAim */) {
            gibbedCastles = [];
            castleGibStates = [];
            loadLevel(state.level);
        }
        modeT = 0;
    }
    // World
    G.useProgram(mainShader);
    G.uniformMatrix4fv(G.getUniformLocation(mainShader, 'a3'), false, m4Ident);
    G.uniformMatrix4fv(G.getUniformLocation(mainShader, 'a2'), false, vp);
    G.uniform3fv(G.getUniformLocation(mainShader, 'ag'), state.pos);
    bindTextureUniforms(mainShader);
    modelGeoDraw(worldGetGeo(), mainShader);
    // Skybox
    G.disable(CULL_FACE);
    G.useProgram(skyShader);
    G.uniformMatrix4fv(G.getUniformLocation(skyShader, 'aa'), false, m4Mul(projectionMat, lookMat));
    modelGeoDraw(skyboxGeo, skyShader);
    G.enable(CULL_FACE);
    // Aim line
    if ((state.mode == 1 /* FirstAim */ || state.mode == 2 /* LaterAim */) && drawCannon) {
        let [predicted, pos] = predictShot(lazyYaw, lazyPitch, state.pos, state.castlesHit);
        ballPos = pos;
        G.bindBuffer(ARRAY_BUFFER, drawRayVertex);
        G.bufferSubData(ARRAY_BUFFER, 0, predicted);
        G.useProgram(aimRayShader);
        G.uniformMatrix4fv(G.getUniformLocation(aimRayShader, 'a3'), false, m4Ident);
        G.uniformMatrix4fv(G.getUniformLocation(aimRayShader, 'a2'), false, vp);
        G.bindBuffer(ARRAY_BUFFER, drawRayVertex);
        let pposLoc = G.getAttribLocation(aimRayShader, 'a1');
        G.enableVertexAttribArray(pposLoc);
        G.vertexAttribPointer(pposLoc, 3, FLOAT, false, 0, 0);
        G.bindBuffer(ELEMENT_ARRAY_BUFFER, drawRayIndex);
        G.drawElements(LINES, 2 * 150, UNSIGNED_SHORT, 0);
        drawBall(vp, aimRayShader, [0, 1, .1]);
    }
    // UI draw
    if (state.mode == 0 /* Menu */ || state.mode == 1 /* FirstAim */) {
        fade = Math.min((modeT / 1000) ** 2, 1);
    }
    else if (state.mode == 4 /* Dead */ || state.mode == 5 /* Win */) {
        fade = 1 - Math.min((modeT / 1000) ** 2, 1);
    }
    else {
        fade = 1;
    }
    let newhash = `${state.castlesHit.length}.${state.ammo}.${state.mode}`;
    if (newhash != uihash) {
        ctx.clearRect(0, 0, C2.width, C2.height);
        ctx.strokeStyle = '#e56b70';
        ctx.fillStyle = '#ffdb63';
        ctx.textAlign = 'center';
        if (state.mode == 0 /* Menu */) {
            ctx.font = 'bold 64px sans-serif';
            ctx.lineWidth = 3;
            drawText("CANNONBOLF", C2.width / 2, 200);
            ctx.font = 'bold 16px sans-serif';
            ctx.lineWidth = .5;
            drawText("CLICK TO START", C2.width / 2, C2.height - 100);
        }
        else if (state.mode == 4 /* Dead */ || state.mode == 5 /* Win */) {
            ctx.font = 'bold 32px sans-serif';
            ctx.lineWidth = 2;
            drawText(state.mode == 4 /* Dead */ ? "SORRY, TRY AGAIN!" : state.level < lastLevel ? (state.ammo > 0 ? "BIRDIE!" : "PAR!") : "THANKS FOR PLAYING!", C2.width / 2, C2.height / 2);
        }
        else {
            ctx.font = 'bold 32px sans-serif';
            ctx.lineWidth = 2;
            ctx.textAlign = 'right';
            drawText(state.ammo + " ⚫", C2.width - 20, C2.height - 25 - 2 * 40);
            drawText(`${state.castlesHit.length}/${worldGetCastles().length} 🏰`, C2.width - 20, C2.height - 25 - 40);
            drawText(`${state.level + 1}/${lastLevel + 1} ⛳`, C2.width - 20, C2.height - 25);
            ctx.textAlign = 'center';
            ctx.font = 'bold 16px sans-serif';
            ctx.lineWidth = .5;
            drawText(hints[state.level], C2.width / 2, C2.height - 20);
        }
    }
    // UI blit
    G.disable(DEPTH_TEST);
    G.useProgram(blitShader);
    G.activeTexture(TEXTURE0);
    G.bindTexture(TEXTURE_2D, ctxtx);
    if (newhash != uihash) {
        G.texImage2D(TEXTURE_2D, 0, RGBA, RGBA, UNSIGNED_BYTE, C2);
    }
    G.uniform1i(G.getUniformLocation(blitShader, 'a4'), 0);
    G.uniform4f(G.getUniformLocation(blitShader, 'a5'), C2.width, C2.height, fade, ~~(state.mode == 5 /* Win */ && state.level == lastLevel));
    let posLoc = G.getAttribLocation(blitShader, 'a1');
    G.bindBuffer(ARRAY_BUFFER, blitTriBuffer);
    G.enableVertexAttribArray(posLoc);
    G.vertexAttribPointer(posLoc, 2, FLOAT, false, 0, 0);
    G.drawArrays(TRIANGLES, 0, 6);
    G.enable(DEPTH_TEST);
    uihash = newhash;
};
let drawText = (txt, x, y) => {
    ctx.fillText(txt, x, y);
    ctx.strokeText(txt, x, y);
};
let drawBall = (vp, shader, color) => {
    let modelMat = m4Mul(m4Mul(m4Translate(ballPos), m4Scale(0.2)), ballRot);
    G.useProgram(shader);
    G.uniformMatrix4fv(G.getUniformLocation(shader, 'a3'), false, modelMat);
    G.uniformMatrix4fv(G.getUniformLocation(shader, 'a2'), false, vp);
    if (shader == aimRayShader) {
        G.uniform3fv(G.getUniformLocation(shader, 'a0'), color);
    }
    else {
        bindTextureUniforms(shader);
    }
    modelGeoDraw(playerGeo, shader);
};

let prevState = gameStateNew();
let curState = gameStateNew();
[document.body, CC].map((elem) => {
    let style = elem.style;
    style.overflow = 'hidden';
    style.margin = 0;
    style.width = '100%';
    style.height = '100%';
    style.cursor = 'pointer';
});
let accTime = 0;
let prevNow = 0;
let accTickInputs = inputsNew();
let frame = (now) => {
    requestAnimationFrame(frame);
    prevNow = prevNow || now;
    let dt = Math.min(now - prevNow, 1000);
    let frameInputs = inputsConsumeFrame();
    let didRunTick = False;
    accTime += dt;
    prevNow = now;
    inputsAdd(accTickInputs, frameInputs);
    while (accTime > 33) {
        didRunTick = True;
        accTime -= 33;
        prevState = curState;
        curState = gameStateTick(curState, accTickInputs);
        accTickInputs.mouseAccX = accTickInputs.mouseAccY = 0;
    }
    if (didRunTick) {
        accTickInputs = inputsNew();
    }
    renderGame(accTickInputs, gameStateLerp(prevState, curState, accTime / 33), dt);
};
window.onresize = resize;
resize();
loadLevel(START_LEVEL);
frame(0);
