let cameras = [
    {
            id: 0,
        img_name: "00001",
        width: 1959,
        height: 1090,
        position: [2, 0, 2],
        rotation: [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ],
        fy: 1164.6601287484507,
        fx: 1159.5880733038064,
    },
    
];

let camera = cameras[0];

function getProjectionMatrix(fx, fy, width, height) {
    const znear = 0.2;
    const zfar = 200;
    return [
        [(2 * fx) / width, 0, 0, 0],
        [0, -(2 * fy) / height, 0, 0],
        [0, 0, zfar / (zfar - znear), 1],
        [0, 0, -(zfar * znear) / (zfar - znear), 0],
    ].flat();
}

function getViewMatrix(camera) {
    const R = camera.rotation.flat();
    const t = camera.position;
    const camToWorld = [
        [R[0], R[1], R[2], 0],
        [R[3], R[4], R[5], 0],
        [R[6], R[7], R[8], 0],
        [
            -t[0] * R[0] - t[1] * R[3] - t[2] * R[6],
            -t[0] * R[1] - t[1] * R[4] - t[2] * R[7],
            -t[0] * R[2] - t[1] * R[5] - t[2] * R[8],
            1,
        ],
    ].flat();
    return camToWorld;
}
// function translate4(a, x, y, z) {
//     return [
//         ...a.slice(0, 12),
//         a[0] * x + a[4] * y + a[8] * z + a[12],
//         a[1] * x + a[5] * y + a[9] * z + a[13],
//         a[2] * x + a[6] * y + a[10] * z + a[14],
//         a[3] * x + a[7] * y + a[11] * z + a[15],
//     ];
// }

function multiply4(a, b) {
    return [
        b[0] * a[0] + b[1] * a[4] + b[2] * a[8] + b[3] * a[12],
        b[0] * a[1] + b[1] * a[5] + b[2] * a[9] + b[3] * a[13],
        b[0] * a[2] + b[1] * a[6] + b[2] * a[10] + b[3] * a[14],
        b[0] * a[3] + b[1] * a[7] + b[2] * a[11] + b[3] * a[15],
        b[4] * a[0] + b[5] * a[4] + b[6] * a[8] + b[7] * a[12],
        b[4] * a[1] + b[5] * a[5] + b[6] * a[9] + b[7] * a[13],
        b[4] * a[2] + b[5] * a[6] + b[6] * a[10] + b[7] * a[14],
        b[4] * a[3] + b[5] * a[7] + b[6] * a[11] + b[7] * a[15],
        b[8] * a[0] + b[9] * a[4] + b[10] * a[8] + b[11] * a[12],
        b[8] * a[1] + b[9] * a[5] + b[10] * a[9] + b[11] * a[13],
        b[8] * a[2] + b[9] * a[6] + b[10] * a[10] + b[11] * a[14],
        b[8] * a[3] + b[9] * a[7] + b[10] * a[11] + b[11] * a[15],
        b[12] * a[0] + b[13] * a[4] + b[14] * a[8] + b[15] * a[12],
        b[12] * a[1] + b[13] * a[5] + b[14] * a[9] + b[15] * a[13],
        b[12] * a[2] + b[13] * a[6] + b[14] * a[10] + b[15] * a[14],
        b[12] * a[3] + b[13] * a[7] + b[14] * a[11] + b[15] * a[15],
    ];
}

function invert4(a) {
    let b00 = a[0] * a[5] - a[1] * a[4];
    let b01 = a[0] * a[6] - a[2] * a[4];
    let b02 = a[0] * a[7] - a[3] * a[4];
    let b03 = a[1] * a[6] - a[2] * a[5];
    let b04 = a[1] * a[7] - a[3] * a[5];
    let b05 = a[2] * a[7] - a[3] * a[6];
    let b06 = a[8] * a[13] - a[9] * a[12];
    let b07 = a[8] * a[14] - a[10] * a[12];
    let b08 = a[8] * a[15] - a[11] * a[12];
    let b09 = a[9] * a[14] - a[10] * a[13];
    let b10 = a[9] * a[15] - a[11] * a[13];
    let b11 = a[10] * a[15] - a[11] * a[14];
    let det =
        b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
    if (!det) return null;
    return [
        (a[5] * b11 - a[6] * b10 + a[7] * b09) / det,
        (a[2] * b10 - a[1] * b11 - a[3] * b09) / det,
        (a[13] * b05 - a[14] * b04 + a[15] * b03) / det,
        (a[10] * b04 - a[9] * b05 - a[11] * b03) / det,
        (a[6] * b08 - a[4] * b11 - a[7] * b07) / det,
        (a[0] * b11 - a[2] * b08 + a[3] * b07) / det,
        (a[14] * b02 - a[12] * b05 - a[15] * b01) / det,
        (a[8] * b05 - a[10] * b02 + a[11] * b01) / det,
        (a[4] * b10 - a[5] * b08 + a[7] * b06) / det,
        (a[1] * b08 - a[0] * b10 - a[3] * b06) / det,
        (a[12] * b04 - a[13] * b02 + a[15] * b00) / det,
        (a[9] * b02 - a[8] * b04 - a[11] * b00) / det,
        (a[5] * b07 - a[4] * b09 - a[6] * b06) / det,
        (a[0] * b09 - a[1] * b07 + a[2] * b06) / det,
        (a[13] * b01 - a[12] * b03 - a[14] * b00) / det,
        (a[8] * b03 - a[9] * b01 + a[10] * b00) / det,
    ];
}

function rotate4(a, rad, x, y, z) {
    let len = Math.hypot(x, y, z);
    x /= len;
    y /= len;
    z /= len;
    let s = Math.sin(rad);
    let c = Math.cos(rad);
    let t = 1 - c;
    let b00 = x * x * t + c;
    let b01 = y * x * t + z * s;
    let b02 = z * x * t - y * s;
    let b10 = x * y * t - z * s;
    let b11 = y * y * t + c;
    let b12 = z * y * t + x * s;
    let b20 = x * z * t + y * s;
    let b21 = y * z * t - x * s;
    let b22 = z * z * t + c;
    return [
        a[0] * b00 + a[4] * b01 + a[8] * b02,
        a[1] * b00 + a[5] * b01 + a[9] * b02,
        a[2] * b00 + a[6] * b01 + a[10] * b02,
        a[3] * b00 + a[7] * b01 + a[11] * b02,
        a[0] * b10 + a[4] * b11 + a[8] * b12,
        a[1] * b10 + a[5] * b11 + a[9] * b12,
        a[2] * b10 + a[6] * b11 + a[10] * b12,
        a[3] * b10 + a[7] * b11 + a[11] * b12,
        a[0] * b20 + a[4] * b21 + a[8] * b22,
        a[1] * b20 + a[5] * b21 + a[9] * b22,
        a[2] * b20 + a[6] * b21 + a[10] * b22,
        a[3] * b20 + a[7] * b21 + a[11] * b22,
        ...a.slice(12, 16),
    ];
}

function translate4(a, x, y, z) {
    return [
        ...a.slice(0, 12),
        a[0] * x + a[4] * y + a[8] * z + a[12],
        a[1] * x + a[5] * y + a[9] * z + a[13],
        a[2] * x + a[6] * y + a[10] * z + a[14],
        a[3] * x + a[7] * y + a[11] * z + a[15],
    ];
}

function createWorker(self) {
    let buffer;
    let vertexCount = 0;
    let viewProj;
    // 6*4 + 4 + 4 = 8*4
    // XYZ - Position (Float32)
    // XYZ - Scale (Float32)
    // RGBA - colors (uint8)
    // IJKL - quaternion/rot (uint8)
    const rowLength = 3 * 4 + 3 * 4 + 4 + 4;
    let lastProj = [];
    let depthIndex = new Uint32Array();
    let lastVertexCount = 0;

    var _floatView = new Float32Array(1);
    var _int32View = new Int32Array(_floatView.buffer);

    function floatToHalf(float) {
        _floatView[0] = float;
        var f = _int32View[0];

        var sign = (f >> 31) & 0x0001;
        var exp = (f >> 23) & 0x00ff;
        var frac = f & 0x007fffff;

        var newExp;
        if (exp == 0) {
            newExp = 0;
        } else if (exp < 113) {
            newExp = 0;
            frac |= 0x00800000;
            frac = frac >> (113 - exp);
            if (frac & 0x01000000) {
                newExp = 1;
                frac = 0;
            }
        } else if (exp < 142) {
            newExp = exp - 112;
        } else {
            newExp = 31;
            frac = 0;
        }

        return (sign << 15) | (newExp << 10) | (frac >> 13);
    }

    function packHalf2x16(x, y) {
        return (floatToHalf(x) | (floatToHalf(y) << 16)) >>> 0;
    }

    function generateTexture() {
        if (!buffer) return;
        const f_buffer = new Float32Array(buffer);
        const u_buffer = new Uint8Array(buffer);

        var texwidth = 1024 * 2; // Set to your desired width
        var texheight = Math.ceil((2 * vertexCount) / texwidth); // Set to your desired height
        var texdata = new Uint32Array(texwidth * texheight * 4); // 4 components per pixel (RGBA)
        var texdata_c = new Uint8Array(texdata.buffer);
        var texdata_f = new Float32Array(texdata.buffer);

        // Here we convert from a .splat file buffer into a texture
        // With a little bit more foresight perhaps this texture file
        // should have been the native format as it'd be very easy to
        // load it into webgl.
        for (let i = 0; i < vertexCount; i++) {
            // x, y, z
            texdata_f[8 * i + 0] = f_buffer[8 * i + 0];
            texdata_f[8 * i + 1] = f_buffer[8 * i + 1];
            texdata_f[8 * i + 2] = f_buffer[8 * i + 2];

            // r, g, b, a
            texdata_c[4 * (8 * i + 7) + 0] = u_buffer[32 * i + 24 + 0];
            texdata_c[4 * (8 * i + 7) + 1] = u_buffer[32 * i + 24 + 1];
            texdata_c[4 * (8 * i + 7) + 2] = u_buffer[32 * i + 24 + 2];
            texdata_c[4 * (8 * i + 7) + 3] = u_buffer[32 * i + 24 + 3];

            // quaternions
            let scale = [
                f_buffer[8 * i + 3 + 0],
                f_buffer[8 * i + 3 + 1],
                f_buffer[8 * i + 3 + 2],
            ];
            let rot = [
                (u_buffer[32 * i + 28 + 0] - 128) / 128,
                (u_buffer[32 * i + 28 + 1] - 128) / 128,
                (u_buffer[32 * i + 28 + 2] - 128) / 128,
                (u_buffer[32 * i + 28 + 3] - 128) / 128,
            ];

            // Compute the matrix product of S and R (M = S * R)
            const M = [
                1.0 - 2.0 * (rot[2] * rot[2] + rot[3] * rot[3]),
                2.0 * (rot[1] * rot[2] + rot[0] * rot[3]),
                2.0 * (rot[1] * rot[3] - rot[0] * rot[2]),

                2.0 * (rot[1] * rot[2] - rot[0] * rot[3]),
                1.0 - 2.0 * (rot[1] * rot[1] + rot[3] * rot[3]),
                2.0 * (rot[2] * rot[3] + rot[0] * rot[1]),

                2.0 * (rot[1] * rot[3] + rot[0] * rot[2]),
                2.0 * (rot[2] * rot[3] - rot[0] * rot[1]),
                1.0 - 2.0 * (rot[1] * rot[1] + rot[2] * rot[2]),
            ].map((k, i) => k * scale[Math.floor(i / 3)]);

            const sigma = [
                M[0] * M[0] + M[3] * M[3] + M[6] * M[6],
                M[0] * M[1] + M[3] * M[4] + M[6] * M[7],
                M[0] * M[2] + M[3] * M[5] + M[6] * M[8],
                M[1] * M[1] + M[4] * M[4] + M[7] * M[7],
                M[1] * M[2] + M[4] * M[5] + M[7] * M[8],
                M[2] * M[2] + M[5] * M[5] + M[8] * M[8],
            ];

            texdata[8 * i + 4] = packHalf2x16(4 * sigma[0], 4 * sigma[1]);
            texdata[8 * i + 5] = packHalf2x16(4 * sigma[2], 4 * sigma[3]);
            texdata[8 * i + 6] = packHalf2x16(4 * sigma[4], 4 * sigma[5]);
        }

        self.postMessage({ texdata, texwidth, texheight }, [texdata.buffer]);
    }

    function runSort(viewProj) {
        if (!buffer) return;
        const f_buffer = new Float32Array(buffer);
        if (lastVertexCount == vertexCount) {
            let dot =
                lastProj[2] * viewProj[2] +
                lastProj[6] * viewProj[6] +
                lastProj[10] * viewProj[10];
            if (Math.abs(dot - 1) < 0.01) {
                return;
            }
        } else {
            generateTexture();
            lastVertexCount = vertexCount;
        }

        console.time("sort");
        let maxDepth = -Infinity;
        let minDepth = Infinity;
        let sizeList = new Int32Array(vertexCount);
        for (let i = 0; i < vertexCount; i++) {
            let depth =
                ((viewProj[2] * f_buffer[8 * i + 0] +
                    viewProj[6] * f_buffer[8 * i + 1] +
                    viewProj[10] * f_buffer[8 * i + 2]) *
                    4096) |
                0;
            sizeList[i] = depth;
            if (depth > maxDepth) maxDepth = depth;
            if (depth < minDepth) minDepth = depth;
        }

        // This is a 16 bit single-pass counting sort
        let depthInv = (256 * 256) / (maxDepth - minDepth);
        let counts0 = new Uint32Array(256 * 256);
        for (let i = 0; i < vertexCount; i++) {
            sizeList[i] = ((sizeList[i] - minDepth) * depthInv) | 0;
            counts0[sizeList[i]]++;
        }
        let starts0 = new Uint32Array(256 * 256);
        for (let i = 1; i < 256 * 256; i++)
            starts0[i] = starts0[i - 1] + counts0[i - 1];
        depthIndex = new Uint32Array(vertexCount);
        for (let i = 0; i < vertexCount; i++)
            depthIndex[starts0[sizeList[i]]++] = i;

        console.timeEnd("sort");

        lastProj = viewProj;
        self.postMessage({ depthIndex, viewProj, vertexCount }, [
            depthIndex.buffer,
        ]);
    }

    function processPlyBuffer(inputBuffer) {
        const ubuf = new Uint8Array(inputBuffer);
        // 10KB ought to be enough for a header...
        const header = new TextDecoder().decode(ubuf.slice(0, 1024 * 10));
        const header_end = "end_header\n";
        const header_end_index = header.indexOf(header_end);
        if (header_end_index < 0)
            throw new Error("Unable to read .ply file header");
        const vertexCount = parseInt(/element vertex (\d+)\n/.exec(header)[1]);
        console.log("Vertex Count", vertexCount);
        let row_offset = 0,
            offsets = {},
            types = {};
        const TYPE_MAP = {
            double: "getFloat64",
            int: "getInt32",
            uint: "getUint32",
            float: "getFloat32",
            short: "getInt16",
            ushort: "getUint16",
            uchar: "getUint8",
        };
        for (let prop of header
            .slice(0, header_end_index)
            .split("\n")
            .filter((k) => k.startsWith("property "))) {
            const [p, type, name] = prop.split(" ");
            const arrayType = TYPE_MAP[type] || "getInt8";
            types[name] = arrayType;
            offsets[name] = row_offset;
            row_offset += parseInt(arrayType.replace(/[^\d]/g, "")) / 8;
        }
        console.log("Bytes per row", row_offset, types, offsets);

        let dataView = new DataView(
            inputBuffer,
            header_end_index + header_end.length,
        );
        let row = 0;
        const attrs = new Proxy(
            {},
            {
                get(target, prop) {
                    if (!types[prop]) throw new Error(prop + " not found");
                    return dataView[types[prop]](
                        row * row_offset + offsets[prop],
                        true,
                    );
                },
            },
        );

        console.time("calculate importance");
        let sizeList = new Float32Array(vertexCount);
        let sizeIndex = new Uint32Array(vertexCount);
        for (row = 0; row < vertexCount; row++) {
            sizeIndex[row] = row;
            if (!types["scale_0"]) continue;
            const size =
                Math.exp(attrs.scale_0) *
                Math.exp(attrs.scale_1) *
                Math.exp(attrs.scale_2);
            const opacity = 1 / (1 + Math.exp(-attrs.opacity));
            sizeList[row] = size * opacity;
        }
        console.timeEnd("calculate importance");

        console.time("sort");
        sizeIndex.sort((b, a) => sizeList[a] - sizeList[b]);
        console.timeEnd("sort");

        // 6*4 + 4 + 4 = 8*4
        // XYZ - Position (Float32)
        // XYZ - Scale (Float32)
        // RGBA - colors (uint8)
        // IJKL - quaternion/rot (uint8)
        const rowLength = 3 * 4 + 3 * 4 + 4 + 4;
        const buffer = new ArrayBuffer(rowLength * vertexCount);

        console.time("build buffer");
        for (let j = 0; j < vertexCount; j++) {
            row = sizeIndex[j];

            const position = new Float32Array(buffer, j * rowLength, 3);
            const scales = new Float32Array(buffer, j * rowLength + 4 * 3, 3);
            const rgba = new Uint8ClampedArray(
                buffer,
                j * rowLength + 4 * 3 + 4 * 3,
                4,
            );
            const rot = new Uint8ClampedArray(
                buffer,
                j * rowLength + 4 * 3 + 4 * 3 + 4,
                4,
            );

            if (types["scale_0"]) {
                const qlen = Math.sqrt(
                    attrs.rot_0 ** 2 +
                        attrs.rot_1 ** 2 +
                        attrs.rot_2 ** 2 +
                        attrs.rot_3 ** 2,
                );

                rot[0] = (attrs.rot_0 / qlen) * 128 + 128;
                rot[1] = (attrs.rot_1 / qlen) * 128 + 128;
                rot[2] = (attrs.rot_2 / qlen) * 128 + 128;
                rot[3] = (attrs.rot_3 / qlen) * 128 + 128;

                scales[0] = Math.exp(attrs.scale_0);
                scales[1] = Math.exp(attrs.scale_1);
                scales[2] = Math.exp(attrs.scale_2);
            } else {
                scales[0] = 0.01;
                scales[1] = 0.01;
                scales[2] = 0.01;

                rot[0] = 255;
                rot[1] = 0;
                rot[2] = 0;
                rot[3] = 0;
            }

            position[0] = attrs.x;
            position[1] = attrs.y;
            position[2] = attrs.z;

            if (types["f_dc_0"]) {
                const SH_C0 = 0.28209479177387814;
                rgba[0] = (0.5 + SH_C0 * attrs.f_dc_0) * 255;
                rgba[1] = (0.5 + SH_C0 * attrs.f_dc_1) * 255;
                rgba[2] = (0.5 + SH_C0 * attrs.f_dc_2) * 255;
            } else {
                rgba[0] = attrs.red;
                rgba[1] = attrs.green;
                rgba[2] = attrs.blue;
            }
            if (types["opacity"]) {
                rgba[3] = (1 / (1 + Math.exp(-attrs.opacity))) * 255;
            } else {
                rgba[3] = 255;
            }
        }
        console.timeEnd("build buffer");
        return buffer;
    }

    const throttledSort = () => {
        if (!sortRunning) {
            sortRunning = true;
            let lastView = viewProj;
            runSort(lastView);
            setTimeout(() => {
                sortRunning = false;
                if (lastView !== viewProj) {
                    throttledSort();
                }
            }, 0);
        }
    };

    let sortRunning;
    self.onmessage = (e) => {
        if (e.data.ply) {
            vertexCount = 0;
            runSort(viewProj);
            buffer = processPlyBuffer(e.data.ply);
            vertexCount = Math.floor(buffer.byteLength / rowLength);
            postMessage({ buffer: buffer });
        } else if (e.data.buffer) {
            buffer = e.data.buffer;
            vertexCount = e.data.vertexCount;
        } else if (e.data.vertexCount) {
            vertexCount = e.data.vertexCount;
        } else if (e.data.view) {
            viewProj = e.data.view;
            throttledSort();
        }
    };
}

const vertexShaderSource = `
#version 300 es
precision highp float;
precision highp int;

uniform highp usampler2D u_texture;
uniform mat4 projection, view;
uniform vec2 focal;
uniform vec2 viewport;

in vec2 position;
in int index;

out vec4 vColor;
out vec2 vPosition;

void main () {
    uvec4 cen = texelFetch(u_texture, ivec2((uint(index) & 0x3ffu) << 1, uint(index) >> 10), 0);
    vec4 cam = view * vec4(uintBitsToFloat(cen.xyz), 1);
    vec4 pos2d = projection * cam;

    float clip = 1.2 * pos2d.w;
    if (pos2d.z < -clip || pos2d.x < -clip || pos2d.x > clip || pos2d.y < -clip || pos2d.y > clip) {
        gl_Position = vec4(0.0, 0.0, 2.0, 1.0);
        return;
    }

    uvec4 cov = texelFetch(u_texture, ivec2(((uint(index) & 0x3ffu) << 1) | 1u, uint(index) >> 10), 0);
    vec2 u1 = unpackHalf2x16(cov.x), u2 = unpackHalf2x16(cov.y), u3 = unpackHalf2x16(cov.z);
    mat3 Vrk = mat3(u1.x, u1.y, u2.x, u1.y, u2.y, u3.x, u2.x, u3.x, u3.y);

    mat3 J = mat3(
        focal.x / cam.z, 0., -(focal.x * cam.x) / (cam.z * cam.z), 
        0., -focal.y / cam.z, (focal.y * cam.y) / (cam.z * cam.z), 
        0., 0., 0.
    );

    mat3 T = transpose(mat3(view)) * J;
    mat3 cov2d = transpose(T) * Vrk * T;

    float mid = (cov2d[0][0] + cov2d[1][1]) / 2.0;
    float radius = length(vec2((cov2d[0][0] - cov2d[1][1]) / 2.0, cov2d[0][1]));
    float lambda1 = mid + radius, lambda2 = mid - radius;

    if(lambda2 < 0.0) return;
    vec2 diagonalVector = normalize(vec2(cov2d[0][1], lambda1 - cov2d[0][0]));
    vec2 majorAxis = min(sqrt(2.0 * lambda1), 1024.0) * diagonalVector;
    vec2 minorAxis = min(sqrt(2.0 * lambda2), 1024.0) * vec2(diagonalVector.y, -diagonalVector.x);

    vColor = clamp(pos2d.z/pos2d.w+1.0, 0.0, 1.0) * vec4((cov.w) & 0xffu, (cov.w >> 8) & 0xffu, (cov.w >> 16) & 0xffu, (cov.w >> 24) & 0xffu) / 255.0;
    vPosition = position;

    vec2 vCenter = vec2(pos2d) / pos2d.w;
    gl_Position = vec4(
        vCenter 
        + position.x * majorAxis / viewport 
        + position.y * minorAxis / viewport, 0.0, 1.0);

}
`.trim();

const fragmentShaderSource = `
#version 300 es
precision highp float;

in vec4 vColor;
in vec2 vPosition;

out vec4 fragColor;

void main () {
    float A = -dot(vPosition, vPosition);
    if (A < -4.0) discard;
    float B = exp(A) * vColor.a;
    fragColor = vec4(B * vColor.rgb, B);
}

`.trim();

async function mergeStreams(reader1, reader2) {
    const { readable, writable } = new TransformStream();

    const writer = writable.getWriter();

    async function pump(reader) {
        while (true) {
            const { done, value } = await reader.read();
            if (done) break;
            writer.write(value);
        }
    }

    // Start pumping both readers
    await pump(reader1);
    await pump(reader2);

    writer.close();

    return readable;
}


let defaultViewMatrix = [
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    -2, 0, -2, 1
];
let viewMatrix = defaultViewMatrix;
let actualViewMatrix = defaultViewMatrix;


function getReaderForFile(filePath) {
    return fetch(filePath)
        .then(response => {
            if (!response.ok) {
                throw new Error('Network response was not ok ' + response.statusText);
            }
            return response.blob();
        })
        .then(blob => {
            return new Promise((resolve, reject) => {
                const reader = new FileReader();
                reader.onload = function(event) {
                    const arrayBuffer = event.target.result;
                    resolve({
                        result: arrayBuffer,
                        read: function () {
                            let position = 0;
                            return {
                                read() {
                                    if (position >= arrayBuffer.byteLength) {
                                        return { done: true };
                                    }
                                    const chunkSize = 1024 * 1024; // 1MB per chunk
                                    const chunk = arrayBuffer.slice(position, position + chunkSize);
                                    position += chunkSize;
                                    return { done: false, value: new Uint8Array(chunk) };
                                }
                            }
                        }
                    });
                };
                reader.onerror = function(error) {
                    reject(error);
                };
                reader.readAsArrayBuffer(blob);
            });
        })
        .catch(error => {
            console.error('Error fetching the local file:', error);
        });
}





async function main() {
    let carousel = false;
    const params = new URLSearchParams(location.search);
    try {
        viewMatrix = JSON.parse(decodeURIComponent(location.hash.slice(1)));
        carousel = false;
    } catch (err) {}

    const url1 = new URL(
        // "nike.splat",
        // location.href,
        params.get("url") || "chunk_0_0.splat",
        "https://ece4500jsplatdata.blob.core.windows.net/wenbobuilding/",
    );
    const req1 = await fetch(url1, {
        headers: {
            'Cache-Control': 'no-store, no-cache, must-revalidate, max-age=0',
            'Pragma': 'no-cache',
            'Expires': '0'
        }
    });

    const url2 = new URL(
        // "nike.splat",
        // location.href,
        params.get("url") || "chunk_-4_0.splat",
        "https://ece4500jsplatdata.blob.core.windows.net/wenbobuilding/",
    );
    const req2 = await fetch(url2, {
        headers: {
            'Cache-Control': 'no-store, no-cache, must-revalidate, max-age=0',
            'Pragma': 'no-cache',
            'Expires': '0'
        }
    });

    const url3 = new URL(
        // "nike.splat",
        // location.href,
        params.get("url") || "chunk_0_-4.splat",
        "https://ece4500jsplatdata.blob.core.windows.net/wenbobuilding/",
    );
    const req3 = await fetch(url3, {
        headers: {
            'Cache-Control': 'no-store, no-cache, must-revalidate, max-age=0',
            'Pragma': 'no-cache',
            'Expires': '0'
        }
    });

    const url4 = new URL(
        // "nike.splat",
        // location.href,
        params.get("url") || "chunk_-4_-4.splat",
        "https://ece4500jsplatdata.blob.core.windows.net/wenbobuilding/",
    );
    const req4 = await fetch(url4, {
        headers: {
            'Cache-Control': 'no-store, no-cache, must-revalidate, max-age=0',
            'Pragma': 'no-cache',
            'Expires': '0'
        }
    });

    const url5 = new URL(
        // "nike.splat",
        // location.href,
        params.get("url") || "chunk_4_0.splat",
        "https://ece4500jsplatdata.blob.core.windows.net/wenbobuilding/",
    );
    const req5 = await fetch(url5, {
        headers: {
            'Cache-Control': 'no-store, no-cache, must-revalidate, max-age=0',
            'Pragma': 'no-cache',
            'Expires': '0'
        }
    });

    const url6 = new URL(
        // "nike.splat",
        // location.href,
        params.get("url") || "chunk_0_4.splat",
        "https://ece4500jsplatdata.blob.core.windows.net/wenbobuilding/",
    );
    const req6 = await fetch(url6, {
        headers: {
            'Cache-Control': 'no-store, no-cache, must-revalidate, max-age=0',
            'Pragma': 'no-cache',
            'Expires': '0'
        }
    });

    const url7 = new URL(
        // "nike.splat",
        // location.href,
        params.get("url") || "chunk_4_4.splat",
        "https://ece4500jsplatdata.blob.core.windows.net/wenbobuilding/",
    );
    const req7 = await fetch(url7, {
        headers: {
            'Cache-Control': 'no-store, no-cache, must-revalidate, max-age=0',
            'Pragma': 'no-cache',
            'Expires': '0'
        }
    });

    const url8 = new URL(
        // "nike.splat",
        // location.href,
        params.get("url") || "chunk_-4_4.splat",
        "https://ece4500jsplatdata.blob.core.windows.net/wenbobuilding/",
    );
    const req8 = await fetch(url8, {
        headers: {
            'Cache-Control': 'no-store, no-cache, must-revalidate, max-age=0',
            'Pragma': 'no-cache',
            'Expires': '0'
        }
    });

    const url9 = new URL(
        // "nike.splat",
        // location.href,
        params.get("url") || "chunk_4_-4.splat",
        "https://ece4500jsplatdata.blob.core.windows.net/wenbobuilding/",
    );
    const req9 = await fetch(url9, {
        headers: {
            'Cache-Control': 'no-store, no-cache, must-revalidate, max-age=0',
            'Pragma': 'no-cache',
            'Expires': '0'
        }
    });

    
    console.log(req1);
    console.log(req2);
    console.log(req3);
    console.log(req4);
    console.log(req5);
    console.log(req6);
    console.log(req7);
    console.log(req8);
    console.log(req9);
    

    if (req1.status != 200)
        throw new Error(req1.status + " Unable to load " + req1.url);
    if (req2.status != 200)
        throw new Error(req2.status + " Unable to load " + req2.url);

    if (!req1.ok) throw new Error(`Failed to fetch ${url1}: ${req1.statusText}`);
    if (!req2.ok) throw new Error(`Failed to fetch ${url2}: ${req2.statusText}`);

    const rowLength = 3 * 4 + 3 * 4 + 4 + 4;

    const reader1 = req1.body.getReader();
    const reader2 = req2.body.getReader();
    const reader3 = req3.body.getReader();
    const reader4 = req4.body.getReader();
    const reader5 = req5.body.getReader();
    const reader6 = req6.body.getReader();
    const reader7 = req7.body.getReader();
    const reader8 = req8.body.getReader();
    const reader9 = req9.body.getReader();
    

    
    let length1 = parseInt(req1.headers.get("content-length"), 10);
    let length2 = parseInt(req2.headers.get("content-length"), 10);
    let length3 = parseInt(req3.headers.get("content-length"), 10);
    let length4 = parseInt(req4.headers.get("content-length"), 10);
    let length5 = parseInt(req5.headers.get("content-length"), 10);
    let length6 = parseInt(req6.headers.get("content-length"), 10);
    let length7 = parseInt(req7.headers.get("content-length"), 10);
    let length8 = parseInt(req8.headers.get("content-length"), 10);
    let length9 = parseInt(req9.headers.get("content-length"), 10);


    let totalLength = (length1 + length2 + length3 + length4 + length5 + length6 + length7 + length8 + length9) * 1;
    let splatData = new Uint8Array(totalLength);

    //let splatData = new Uint8Array(req1.headers.get("content-length") + req2.headers.get("content-length"));

    const downsample =
        splatData.length / rowLength > 500000 ? 1 : 1 / devicePixelRatio;
    console.log(splatData.length / rowLength, downsample);

    const worker = new Worker(
        URL.createObjectURL(
            new Blob(["(", createWorker.toString(), ")(self)"], {
                type: "application/javascript",
            }),
        ),
    );

    const canvas = document.getElementById("canvas");
    const fps = document.getElementById("fps");
    const camid = document.getElementById("camid");

    let projectionMatrix;

    const gl = canvas.getContext("webgl2", {
        antialias: false,
    });

    const vertexShader = gl.createShader(gl.VERTEX_SHADER);
    gl.shaderSource(vertexShader, vertexShaderSource);
    gl.compileShader(vertexShader);
    if (!gl.getShaderParameter(vertexShader, gl.COMPILE_STATUS))
        console.error(gl.getShaderInfoLog(vertexShader));

    const fragmentShader = gl.createShader(gl.FRAGMENT_SHADER);
    gl.shaderSource(fragmentShader, fragmentShaderSource);
    gl.compileShader(fragmentShader);
    if (!gl.getShaderParameter(fragmentShader, gl.COMPILE_STATUS))
        console.error(gl.getShaderInfoLog(fragmentShader));

    const program = gl.createProgram();
    gl.attachShader(program, vertexShader);
    gl.attachShader(program, fragmentShader);
    gl.linkProgram(program);
    gl.useProgram(program);

    if (!gl.getProgramParameter(program, gl.LINK_STATUS))
        console.error(gl.getProgramInfoLog(program));

    gl.disable(gl.DEPTH_TEST); // Disable depth testing

    // Enable blending
    gl.enable(gl.BLEND);
    gl.blendFuncSeparate(
        gl.ONE_MINUS_DST_ALPHA,
        gl.ONE,
        gl.ONE_MINUS_DST_ALPHA,
        gl.ONE,
    );
    gl.blendEquationSeparate(gl.FUNC_ADD, gl.FUNC_ADD);

    const u_projection = gl.getUniformLocation(program, "projection");
    const u_viewport = gl.getUniformLocation(program, "viewport");
    const u_focal = gl.getUniformLocation(program, "focal");
    const u_view = gl.getUniformLocation(program, "view");

    // positions
    const triangleVertices = new Float32Array([-2, -2, 2, -2, 2, 2, -2, 2]);
    const vertexBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, triangleVertices, gl.STATIC_DRAW);
    const a_position = gl.getAttribLocation(program, "position");
    gl.enableVertexAttribArray(a_position);
    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
    gl.vertexAttribPointer(a_position, 2, gl.FLOAT, false, 0, 0);

    var texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, texture);

    var u_textureLocation = gl.getUniformLocation(program, "u_texture");
    gl.uniform1i(u_textureLocation, 0);

    const indexBuffer = gl.createBuffer();
    const a_index = gl.getAttribLocation(program, "index");
    gl.enableVertexAttribArray(a_index);
    gl.bindBuffer(gl.ARRAY_BUFFER, indexBuffer);
    gl.vertexAttribIPointer(a_index, 1, gl.INT, false, 0, 0);
    gl.vertexAttribDivisor(a_index, 1);

    const resize = () => {
        gl.uniform2fv(u_focal, new Float32Array([camera.fx, camera.fy]));

        projectionMatrix = getProjectionMatrix(
            camera.fx,
            camera.fy,
            innerWidth,
            innerHeight,
        );

        gl.uniform2fv(u_viewport, new Float32Array([innerWidth, innerHeight]));

        gl.canvas.width = Math.round(innerWidth / downsample);
        gl.canvas.height = Math.round(innerHeight / downsample);
        gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);

        gl.uniformMatrix4fv(u_projection, false, projectionMatrix);
    };

    window.addEventListener("resize", resize);
    resize();

    worker.onmessage = (e) => {
        if (e.data.buffer) {
            splatData = new Uint8Array(e.data.buffer);
            const blob = new Blob([splatData.buffer], {
                type: "application/octet-stream",
            });
            const link = document.createElement("a");
            link.download = "model.splat";
            link.href = URL.createObjectURL(blob);
            document.body.appendChild(link);
            link.click();
        } else if (e.data.texdata) {
            const { texdata, texwidth, texheight } = e.data;
            // console.log(texdata)
            gl.bindTexture(gl.TEXTURE_2D, texture);
            gl.texParameteri(
                gl.TEXTURE_2D,
                gl.TEXTURE_WRAP_S,
                gl.CLAMP_TO_EDGE,
            );
            gl.texParameteri(
                gl.TEXTURE_2D,
                gl.TEXTURE_WRAP_T,
                gl.CLAMP_TO_EDGE,
            );
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);

            gl.texImage2D(
                gl.TEXTURE_2D,
                0,
                gl.RGBA32UI,
                texwidth,
                texheight,
                0,
                gl.RGBA_INTEGER,
                gl.UNSIGNED_INT,
                texdata,
            );
            gl.activeTexture(gl.TEXTURE0);
            gl.bindTexture(gl.TEXTURE_2D, texture);
        } else if (e.data.depthIndex) {
            const { depthIndex, viewProj } = e.data;
            gl.bindBuffer(gl.ARRAY_BUFFER, indexBuffer);
            gl.bufferData(gl.ARRAY_BUFFER, depthIndex, gl.DYNAMIC_DRAW);
            vertexCount = e.data.vertexCount;
        }
    };

    let activeKeys = [];
	let currentCameraIndex = 0;

    const ensureXZPlaneAlignment = (viewMatrix) => {
        // Extract the current position and direction from the view matrix
        let inv = invert4(viewMatrix);
        let pos = [inv[12], inv[13], inv[14]];  // Camera position
        let dir = [inv[8], inv[9], inv[10]];    // Camera direction
    
        // Force the direction to be parallel to the XZ plane
        dir[1] = 0; // Remove any Y component
        let len = Math.hypot(dir[0], dir[2]);
        dir[0] /= len;
        dir[2] /= len;
    
        // Recompute the up vector to ensure orthogonality
        let up = [0, 1, 0];
        let right = [
            dir[2] * up[1] - dir[1] * up[2],
            dir[0] * up[2] - dir[2] * up[0],
            dir[1] * up[0] - dir[0] * up[1]
        ];
    
        // Reconstruct the view matrix
        inv[0] = right[0]; inv[4] = up[0]; inv[8] = dir[0];
        inv[1] = right[1]; inv[5] = up[1]; inv[9] = dir[1];
        inv[2] = right[2]; inv[6] = up[2]; inv[10] = dir[2];
        inv[12] = pos[0]; inv[13] = pos[1]; inv[14] = pos[2];
    
        return invert4(inv);
    };


    window.addEventListener("keydown", (e) => {
        // if (document.activeElement != document.body) return;
        carousel = false;
        if (!activeKeys.includes(e.code)) activeKeys.push(e.code);
        if (/\d/.test(e.key)) {
            currentCameraIndex = parseInt(e.key)
            camera = camera;
            viewMatrix = getViewMatrix(camera);
        }
		if (['-', '_'].includes(e.key)){
			//currentCameraIndex = (currentCameraIndex + cameras.length - 1) % cameras.length;
			viewMatrix = getViewMatrix(camera);
		}
		if (['+', '='].includes(e.key)){
			//currentCameraIndex = (currentCameraIndex + 1) % cameras.length;
			viewMatrix = getViewMatrix(camera);
		}
        camid.innerText = "cam  " + currentCameraIndex;
        if (e.code == "KeyV") {
            location.hash =
                "#" +
                JSON.stringify(
                    viewMatrix.map((k) => Math.round(k * 100) / 100),
                );
                camid.innerText =""
        } else if (e.code === "KeyP") {
            carousel = false;
            camid.innerText =""
        }
    });
    window.addEventListener("keyup", (e) => {
        activeKeys = activeKeys.filter((k) => k !== e.code);
    });
    window.addEventListener("blur", () => {
        activeKeys = [];
    });

    

    let startX, startY, down;
    canvas.addEventListener("mousedown", (e) => {
        carousel = false;
        e.preventDefault();
        startX = e.clientX;
        startY = e.clientY;
        down = 1
    });
    canvas.addEventListener("contextmenu", (e) => {
        carousel = false;
        e.preventDefault();
        startX = e.clientX;
        startY = e.clientY;
        down = 1;
    });
    
    let mousemove_y = 0;
    let totalVerticalDistance = 0; // Variable to store total vertical distance
    let lastY = null;
    canvas.addEventListener("mousemove", (e) => {
        e.preventDefault();
        if (down == 1) {
            let inv = invert4(viewMatrix);
            let inv2 = invert4(viewMatrix);
            let dx = (5 * (e.clientX - startX)) / innerWidth;
            mousemove_y = (5 * (e.clientY - startY)) / innerHeight;
            let d = 4;

            //inv = translate4(inv, 0, 0, d);
            inv = rotate4(inv, dx, 0, 1, 0);
            //inv = rotate4(inv, -dy, 1, 0, 0);
            //inv = translate4(inv, 0, 0, -d);
            // let postAngle = Math.atan2(inv[0], inv[10])
            // inv = rotate4(inv, postAngle - preAngle, 0, 0, 1)
            // console.log(postAngle)
            viewMatrix = invert4(inv);
            //viewMatrix = ensureXZPlaneAlignment(invert4(inv));
            if (e.buttons === 1) { // Check if the left mouse button is pressed
                if (lastY !== null) {
                    const deltaY = (5 * (e.clientY - startY)) / innerHeight;
                    totalVerticalDistance += deltaY;
                }
                lastY = e.clientY;
            } else {
                lastY = null; // Reset lastY when the mouse button is not pressed
            }
            startX = e.clientX;
            startY = e.clientY;
        } 
    });
    canvas.addEventListener("mouseup", (e) => {
        e.preventDefault();
        down = false;
        startX = 0;
        startY = 0;
    });

    let altX = 0,
        altY = 0;
    canvas.addEventListener(
        "touchstart",
        (e) => {
            e.preventDefault();
            if (e.touches.length === 1) {
                carousel = false;
                startX = e.touches[0].clientX;
                startY = e.touches[0].clientY;
                down = 1;
            } else if (e.touches.length === 2) {
                // console.log('beep')
                carousel = false;
                startX = e.touches[0].clientX;
                altX = e.touches[1].clientX;
                startY = e.touches[0].clientY;
                altY = e.touches[1].clientY;
                down = 1;
            }
        },
        { passive: false },
    );
    canvas.addEventListener(
        "touchmove",
        (e) => {
            e.preventDefault();
            if (e.touches.length === 1 && down) {
                let inv = invert4(viewMatrix);
                let dx = (4 * (e.touches[0].clientX - startX)) / innerWidth;
                let dy = (4 * (e.touches[0].clientY - startY)) / innerHeight;

                let d = 4;
                inv = translate4(inv, 0, 0, d);
                // inv = translate4(inv,  -x, -y, -z);
                // inv = translate4(inv,  x, y, z);
                inv = rotate4(inv, dx, 0, 1, 0);
                inv = rotate4(inv, -dy, 1, 0, 0);
                inv = translate4(inv, 0, 0, -d);

                viewMatrix = invert4(inv);
                actualViewMatrix = invert4(inv);
                startX = e.touches[0].clientX;
                startY = e.touches[0].clientY;
            } else if (e.touches.length === 2) {
                // alert('beep')
                const dtheta =
                    Math.atan2(startY - altY, startX - altX) -
                    Math.atan2(
                        e.touches[0].clientY - e.touches[1].clientY,
                        e.touches[0].clientX - e.touches[1].clientX,
                    );
                const dscale =
                    Math.hypot(startX - altX, startY - altY) /
                    Math.hypot(
                        e.touches[0].clientX - e.touches[1].clientX,
                        e.touches[0].clientY - e.touches[1].clientY,
                    );
                const dx =
                    (e.touches[0].clientX +
                        e.touches[1].clientX -
                        (startX + altX)) /
                    2;
                const dy =
                    (e.touches[0].clientY +
                        e.touches[1].clientY -
                        (startY + altY)) /
                    2;
                let inv = invert4(viewMatrix);
                // inv = translate4(inv,  0, 0, d);
                inv = rotate4(inv, dtheta, 0, 0, 1);

                inv = translate4(inv, -dx / innerWidth, -dy / innerHeight, 0);

                // let preY = inv[13];
                inv = translate4(inv, 0, 0, 3 * (1 - dscale));
                // inv[13] = preY;

                viewMatrix = invert4(inv);

                startX = e.touches[0].clientX;
                altX = e.touches[1].clientX;
                startY = e.touches[0].clientY;
                altY = e.touches[1].clientY;
            }
        },
        { passive: false },
    );
    canvas.addEventListener(
        "touchend",
        (e) => {
            e.preventDefault();
            down = false;
            startX = 0;
            startY = 0;
        },
        { passive: false },
    );

    let jumpDelta = 0;
    let vertexCount = 0;

    let lastFrame = 0;
    let avgFps = 0;
    let start = 0;

    

    let leftGamepadTrigger, rightGamepadTrigger;

    
    //function updateCoordinatesDisplay(position) {
    //    const coordinatesElement = document.getElementById('coordinates');
    //    coordinatesElement.innerText = `Coordinates: (${position[0].toFixed(2)}, ${position[1].toFixed(2)}, ${position[2].toFixed(2)})`;
    //}
    
    function getCurrentPosition(viewMatrix) {
        // Extract position from the view matrix
        let inv = invert4(viewMatrix);
        return [inv[12], inv[13], inv[14]]; // x, y, z coordinates
    }

    function updateCoordinatesDisplay(position, rotation, chunkX, chunkZ) {
        const coordinatesElement = document.getElementById('coordinates');
        coordinatesElement.innerText = `Coordinates: (${position[0].toFixed(2)}, ${position[1].toFixed(2)}, ${position[2].toFixed(2)}), Rotation: (${rotation[0].toFixed(2)}, ${rotation[1].toFixed(2)}, ${rotation[2].toFixed(2)}), Chunk Coordinates: (${chunkX.toFixed(0)}, ${chunkZ.toFixed(0)})`;
    }
    
    function getCurrentPositionAndRotation(viewMatrix) {
        // Extract position from the view matrix
        let inv = invert4(viewMatrix);
        let position = [inv[12], inv[13], inv[14]]; // x, y, z coordinates
    
        // Extract rotation angles from the view matrix
        let rotation = [
            Math.atan2(inv[6], inv[10]), // Rotation around x-axis
            Math.asin(-inv[2]),          // Rotation around y-axis
            Math.atan2(inv[1], inv[0])   // Rotation around z-axis
        ];
    
        // Convert rotation angles from radians to degrees
        rotation = rotation.map(angle => angle * (180 / Math.PI));
    
        return { position, rotation };
    }
    
    function getChunkCoordinates(viewMatrix) {
        // Extract the current position from the view matrix
        let inv = invert4(viewMatrix);
        let position = [inv[12], inv[13], inv[14]]; // Camera position
    
        // Assuming chunk size is 1 (adjust this if your chunks are larger)
        let chunkSize = 4;
        
        // Calculate the lower-left corner of the current chunk
        let chunkX = Math.floor(position[0] / chunkSize) * chunkSize;
        
        let chunkZ = Math.floor(position[2] / chunkSize) * chunkSize;
    
        return {chunkX, chunkZ};
    }


    let lastFrameTime = 0;
    const fpsInterval = 1000 / 30;

    let heightCrossed = false; // Flag to track if height has crossed 5
    let lastHeight = null;
    let lastchunkX = null;
    let lastchunkZ = null;
    let isFiltering = false;

    
    let buttonStates = {
        forward: false,
        backward: false,
        left: false,
        right: false,
        up: false,
        down: false
    };
    function moveCamera(direction, state) {
        buttonStates[direction] = state;
    }

    // Event listeners for the buttons
    document.getElementById('moveForward').addEventListener('mousedown', () => moveCamera('forward', true));
    document.getElementById('moveForward').addEventListener('mouseup', () => moveCamera('forward', false));
    document.getElementById('moveBackward').addEventListener('mousedown', () => moveCamera('backward', true));
    document.getElementById('moveBackward').addEventListener('mouseup', () => moveCamera('backward', false));
    document.getElementById('moveLeft').addEventListener('mousedown', () => moveCamera('left', true));
    document.getElementById('moveLeft').addEventListener('mouseup', () => moveCamera('left', false));
    document.getElementById('moveRight').addEventListener('mousedown', () => moveCamera('right', true));
    document.getElementById('moveRight').addEventListener('mouseup', () => moveCamera('right', false));
    document.getElementById('moveUp').addEventListener('mousedown', () => moveCamera('up', true));
    document.getElementById('moveUp').addEventListener('mouseup', () => moveCamera('up', false));
    document.getElementById('moveDown').addEventListener('mousedown', () => moveCamera('down', true));
    document.getElementById('moveDown').addEventListener('mouseup', () => moveCamera('down', false));


    const frame = async (now) => {
        let elapsed = now - lastFrameTime;

        if (elapsed >= fpsInterval) {
            lastFrameTime = now - (elapsed % fpsInterval);
            let inv = invert4(viewMatrix);
            let shiftKey = activeKeys.includes("Shift") || activeKeys.includes("ShiftLeft") || activeKeys.includes("ShiftRight")

            
            
            
            if (activeKeys.includes("ArrowUp") | buttonStates.up) {
                inv = translate4(inv, 0, -0.03, 0);  // 向上移动
            }
            if (activeKeys.includes("ArrowDown") | buttonStates.down) {
                inv = translate4(inv, 0, 0.03, 0);  // 向下移动
            }
            if (activeKeys.includes("KeyW") | buttonStates.forward) {
                inv = translate4(inv, 0, 0, 0.025);  // 向前移动
            }
            if (activeKeys.includes("KeyS") | buttonStates.backward) {
                inv = translate4(inv, 0, 0, -0.025);  // 向后移动
            }
            if (activeKeys.includes("KeyA") | buttonStates.left) {
                inv = translate4(inv, -0.015, 0, 0);  // 向左移动
            }
            if (activeKeys.includes("KeyD") | buttonStates.right) {
                inv = translate4(inv, 0.015, 0, 0);  // 向右移动
            }
            if (activeKeys.includes("KeyQ")) inv = rotate4(inv, 0.01, 0, 0, 1);
            if (activeKeys.includes("KeyE")) inv = rotate4(inv, -0.01, 0, 0, 1);

            
            let isJumping = activeKeys.includes("Space");
            //let isMoving  = 

            
            //viewMatrix = ensureXZPlaneAlignment(invert4(inv));
            viewMatrix = invert4(inv);

            let { position, rotation } = getCurrentPositionAndRotation(viewMatrix);
            

            let {chunkX,chunkZ}  = getChunkCoordinates(viewMatrix);
            updateCoordinatesDisplay(position, rotation, chunkX, chunkZ);

            if (isDataLoaded && lastHeight !== null && lastHeight <= 5 && position[1] > 5 && !heightCrossed) {
                filterData();
                heightCrossed = true; // Set flag to true to prevent further filtering
            }
            
            

            if (chunkX !== lastchunkX || lastchunkZ !== chunkZ){
                let furchunks = getfurChunks(chunkX, chunkZ, lastchunkX, lastchunkZ);
                
                if(lastchunkX !== null){  
                    filterchunkdata(chunkX, chunkZ, lastchunkX, lastchunkZ);
                }
                

                let newchunks = getnewChunks(chunkX, chunkZ, lastchunkX, lastchunkZ);
                for (const currentChunkCoordinates of newchunks){
                    try {
                        const newChunkData = await loadChunkData(currentChunkCoordinates.x, currentChunkCoordinates.z);
                        
                        // Merge new data with existing data
                        const combinedDataBuffer =  mergeArrayBuffers(splatData.buffer, newChunkData);
                        const combinedData = new Uint8Array(combinedDataBuffer);
                        worker.postMessage({ buffer: combinedData.buffer, vertexCount: Math.floor(combinedData.length / rowLength) });
                        console.log(`New chunk loaded successfully at (${currentChunkCoordinates.x}, ${currentChunkCoordinates.z})`);
                        splatData = combinedData;
                    } catch (error) {
                        console.error('Failed to load new chunk data:', error);
                    }
                }
                
            }

            
            lastHeight = position[1];
            lastchunkX = chunkX;
            lastchunkZ = chunkZ;

            if(heightCrossed && position[1] < 5){
                heightCrossed = false;
            }

            if (carousel) {
                let inv = invert4(defaultViewMatrix);

                const t = Math.sin((Date.now() - start) / 5000);
                inv = translate4(inv, 2.5 * t, 0, 6 * (1 - Math.cos(t)));
                inv = rotate4(inv, -0.6 * t, 0, 1, 0);

                viewMatrix = invert4(inv);
            }

            let move_y_delta = 0 

            

            let inv2 = invert4(viewMatrix);
            //inv2 = translate4(inv2, 0, -jumpDelta, 0);
            inv2 = rotate4(inv2, -0.5*totalVerticalDistance, 1, 0, 0);
            actualViewMatrix = invert4(inv2);

            const viewProj = multiply4(projectionMatrix, actualViewMatrix);
            worker.postMessage({ view: viewProj });

            const currentFps = 1000 / (now - lastFrame) || 0;
            avgFps = avgFps * 0.9 + currentFps * 0.1;

            if (vertexCount > 0) {
                document.getElementById("spinner").style.display = "none";
                gl.uniformMatrix4fv(u_view, false, actualViewMatrix);
                gl.clear(gl.COLOR_BUFFER_BIT);
                gl.drawArraysInstanced(gl.TRIANGLE_FAN, 0, 4, vertexCount);
            } else {
                gl.clear(gl.COLOR_BUFFER_BIT);
                document.getElementById("spinner").style.display = "";
                start = Date.now() + 2000;
            }
            const progress = (100 * vertexCount) / (splatData.length / rowLength);
            if (progress < 100) {
                document.getElementById("progress").style.width = progress + "%";
            } else {
                document.getElementById("progress").style.display = "none";
            }
            fps.innerText = Math.round(avgFps) + " fps";
            if (isNaN(currentCameraIndex)){
                camid.innerText = "";
            }
            lastFrame = now;
            requestAnimationFrame(frame);
        } else {
            requestAnimationFrame(frame);
        }
    };

    frame();


    function getnewChunks(chunkX, chunkZ, lastchunkX, lastchunkZ) {
        let surroundingChunks = [];
        for (let dx = -1; dx <= 1; dx++) {
            for (let dz = -1; dz <= 1; dz++) {
                if(chunkX + dx*4 < lastchunkX - 4 || chunkX + dx*4 > lastchunkX + 4 || chunkZ + dz*4 < lastchunkZ - 4 || chunkZ + dz*4 > lastchunkZ + 4){
                    surroundingChunks.push({ x: chunkX + 4*dx, z: chunkZ + 4*dz });
                }
            }
        }
        
        return surroundingChunks;
    }

    function getfurChunks(chunkX, chunkZ, lastchunkX, lastchunkZ) {
        let furChunks = [];
        for (let dx = -1; dx <= 1; dx++) {
            for (let dz = -1; dz <= 1; dz++) {
                if(chunkX + dx*4 < lastchunkX - 4 || chunkX + dx*4 > lastchunkX + 4 || chunkZ + dz*4 < lastchunkZ - 4 || chunkZ + dz*4 > lastchunkZ + 4){
                    furChunks.push({ x: lastchunkX - dx*4, z: lastchunkZ - dz*4 });
                }
            }
        }
        return furChunks;
    }


    function getChunkUrl(x, z) {
        return `https://ece4500jsplatdata.blob.core.windows.net/wenbobuilding/chunk_${x}_${z}.splat`;
    }

    async function filterchunkdata(x,z,last_x,last_z) {
        // Assuming `splatData` is the array containing your data points
        const rowLength = 3 * 4 + 3 * 4 + 4 + 4;
        let filteredData = [];
        let floatView = new Float32Array(splatData.buffer);
    
        for (let i = 0; i < floatView.length; i += rowLength / 4) {
            if (x == last_x && z > last_z) { // Check if height (y) >= 0
                if (floatView[i+2] > last_z){
                    filteredData.push(...floatView.slice(i, i + rowLength / 4));
                }
            }else if(x == last_x && z < last_z){
                if (floatView[i+2] < last_z+4){
                    filteredData.push(...floatView.slice(i, i + rowLength / 4));
                }
            }else if(x > last_x && z == last_z){
                if (floatView[i] > last_x){
                    filteredData.push(...floatView.slice(i, i + rowLength / 4));
                }
            }else if(x < last_x && z == last_z){
                if (floatView[i] < last_x+4){
                    filteredData.push(...floatView.slice(i, i + rowLength / 4));
                }
            }

           
        }
    
        splatData = new Uint8Array(new ArrayBuffer(filteredData.length * 4));
        new Float32Array(splatData.buffer).set(filteredData);
    
        worker.postMessage({
            buffer: splatData.buffer,
            vertexCount: Math.floor(splatData.length / rowLength),
        });
    }

    async function loadChunkData(chunkX, chunkZ) {
        const chunkUrl = getChunkUrl(chunkX, chunkZ);
        const response = await fetch(chunkUrl);
        if (!response.ok) {
            throw new Error(`Failed to load chunk at ${chunkX}, ${chunkZ}`);
        }
        const reader = response.body.getReader();
        const length = parseInt(response.headers.get("content-length"), 10);
        let chunkData = new Uint8Array(length);
        let bytesRead = 0;
    
        while (true) {
            const { done, value } = await reader.read();
            if (done) break;
            chunkData.set(value, bytesRead);
            bytesRead += value.length;
        }
    
        return chunkData.buffer;
    }

    function mergeArrayBuffers(buffer1, buffer2) {
        const tmp = new Uint8Array(buffer1.byteLength + buffer2.byteLength);
        tmp.set(new Uint8Array(buffer1), 0);
        tmp.set(new Uint8Array(buffer2), buffer1.byteLength);
        return tmp.buffer;
    }


    const selectFile = (file) => {
        const fr = new FileReader();
        if (/\.json$/i.test(file.name)) {
            fr.onload = () => {
                cameras = JSON.parse(fr.result);
                viewMatrix = getViewMatrix(cameras[0]);
                projectionMatrix = getProjectionMatrix(
                    camera.fx / downsample,
                    camera.fy / downsample,
                    canvas.width,
                    canvas.height,
                );
                gl.uniformMatrix4fv(u_projection, false, projectionMatrix);

                console.log("Loaded Cameras");
            };
            fr.readAsText(file);
        } else {
            stopLoading = true;
            fr.onload = () => {
                splatData = new Uint8Array(fr.result);
                console.log("Loaded", Math.floor(splatData.length / rowLength));

                if (
                    splatData[0] == 112 &&
                    splatData[1] == 108 &&
                    splatData[2] == 121 &&
                    splatData[3] == 10
                ) {
                    // ply file magic header means it should be handled differently
                    worker.postMessage({ ply: splatData.buffer });
                } else {
                    worker.postMessage({
                        buffer: splatData.buffer,
                        vertexCount: Math.floor(splatData.length / rowLength),
                    });
                }
            };
            fr.readAsArrayBuffer(file);
        }
    };

    window.addEventListener("hashchange", (e) => {
        try {
            viewMatrix = JSON.parse(decodeURIComponent(location.hash.slice(1)));
            carousel = false;
        } catch (err) {}
    });

    const preventDefault = (e) => {
        e.preventDefault();
        e.stopPropagation();
    };
    document.addEventListener("dragenter", preventDefault);
    document.addEventListener("dragover", preventDefault);
    document.addEventListener("dragleave", preventDefault);
    document.addEventListener("drop", (e) => {
        e.preventDefault();
        e.stopPropagation();
        selectFile(e.dataTransfer.files[0]);
    });

    document.getElementById('uploadButton').addEventListener('click', () => {
        document.getElementById('fileInput').click();
    });
    
    document.getElementById('fileInput').addEventListener('change', (e) => {
        const file = e.target.files[0];
        if (file) {
            selectFile(file);
        }
    });

    

    let isDataLoaded = false;
    let bytesRead = 0;
    let lastVertexCount = -1;
    let stopLoading = false;






    const loadData = async (reader) => {
        while (true) {
            const { done, value } = await reader.read();
            if (done || stopLoading) break;

            splatData.set(value, bytesRead);
            bytesRead += value.length;

            if (vertexCount > lastVertexCount) {
                worker.postMessage({
                    buffer: splatData.buffer,
                    vertexCount: Math.floor(bytesRead / rowLength),
                });
                lastVertexCount = vertexCount;
            }
        }
    };

    const loadAllData = async () => {
        //loadData(reader); 
        //loadData(reader2);
        await Promise.all([loadData(reader1)]);
        await Promise.all([loadData(reader2)]);
        await Promise.all([loadData(reader3)]);
        await Promise.all([loadData(reader4)]);
        await Promise.all([loadData(reader5)]);
        await Promise.all([loadData(reader6)]);
        await Promise.all([loadData(reader7)]);
        await Promise.all([loadData(reader8)]);
        await Promise.all([loadData(reader9)]);
        
        worker.postMessage({
            buffer: splatData.buffer,
            vertexCount: Math.floor(bytesRead / rowLength),
        });
        
        isDataLoaded = true;
    };

    loadAllData();


    function filterData() {
        // Assuming `splatData` is the array containing your data points
        const rowLength = 3 * 4 + 3 * 4 + 4 + 4;
        let filteredData = [];
        let floatView = new Float32Array(splatData.buffer);
    
        for (let i = 0; i < floatView.length; i += rowLength / 4) {
            if (floatView[i + 1] >= 0) { // Check if height (y) >= 0
                filteredData.push(...floatView.slice(i, i + rowLength / 4));
            }

           
        }
    
        splatData = new Uint8Array(new ArrayBuffer(filteredData.length * 4));
        new Float32Array(splatData.buffer).set(filteredData);
    
        worker.postMessage({
            buffer: splatData.buffer,
            vertexCount: Math.floor(splatData.length / rowLength),
        });
    }
    

    if (!stopLoading)
        worker.postMessage({
            buffer: splatData.buffer,
            vertexCount: Math.floor(bytesRead / rowLength),
        });
    
    

}


window.addEventListener('beforeunload', (event) => {
    // Release any resources or caches here
    if (worker) {
        worker.terminate();
    }
    gl.getExtension('WEBGL_lose_context').loseContext();
    console.log("Resources and caches have been released.");

    // Cancel any ongoing requests or processes
    if (ongoingRequests) {
        ongoingRequests.forEach(controller => controller.abort());
    }
    if (processes) {
        processes.forEach(process => process.kill());
    }
    console.log("All ongoing requests and processes have been terminated.");
});


main().catch((err) => {
    document.getElementById("spinner").style.display = "none";
    document.getElementById("message").innerText = err.toString();
});
