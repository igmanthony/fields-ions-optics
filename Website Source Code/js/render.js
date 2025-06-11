import {glob, turbo} from './globals.js';
import m4 from './libraries/m4.js';
import {
    createProgramInfo,
    createBufferInfoFromArrays,
    setUniforms,
    setBuffersAndAttributes,
} from './libraries/twgl.js';

const VERTEX_SHADER = `attribute vec4 position;
attribute vec3 normal;
attribute vec2 texcoord;
uniform mat4 projection;
uniform mat4 modelView;
varying vec3 v_normal;
varying vec2 v_texcoord;
void main() {
  gl_Position = projection * modelView * position;
  v_normal = mat3(modelView) * normal;
  v_texcoord = texcoord;
}
`;

const FRAGMENT_SHADER = `precision highp float;
varying vec3 v_normal;
varying vec2 v_texcoord;
varying float v_modelId;
uniform sampler2D u_texture;
void main() {
  vec3 lightDirection = normalize(vec3(1, 2, -3));  // arbitrary light direction
  float l = dot(lightDirection, normalize(v_normal)) * .5 + .5;
  highp vec4 texelColor = texture2D(u_texture, v_texcoord);
  gl_FragColor = vec4(texelColor.rgb * l, texelColor.a);
}
`;


function degToRad(d) { return d * Math.PI / 180; }

function generateFaceNormals(indices, positions) {
    const faceNormals = [];
    for (let i = 0; i < indices.length / 3; i++) {
        const n1 = indices[i * 3] * 3;
        const n2 = indices[i * 3 + 1] * 3;
        const n3 = indices[i * 3 + 2] * 3;
        const v1 = positions.slice(n1, n1 + 3);
        const v2 = positions.slice(n2, n2 + 3);
        const v3 = positions.slice(n3, n3 + 3);
        faceNormals.push(m4.normalize(m4.cross(m4.subtractVectors(v1, v2), m4.subtractVectors(v3, v2))));
    }
    return faceNormals;
}

function generateVertIndices(positions) {
    let tempVerts = {};
    let tempVertNdx = 0;
    const vertIndices = [];
    for (let i = 0; i < positions.length; i++) {
        const offset = i * 3;
        const vert = positions.slice(offset, offset + 3);
        const vertId = vert[0] + "," + vert[1] + "," + vert[2];
        if (tempVerts[vertId] === undefined) {
            tempVerts[vertId] = tempVertNdx++;
        }
        vertIndices.push(tempVerts[vertId]);
    }
    return vertIndices;
}

function generateVertFaces(vertIndices, indices) {
    const vertFaces = [];
    for (let i = 0; i < indices.length; ++i) {
        const sharedNdx = vertIndices[indices[i]];
        vertFaces[sharedNdx] = (!vertFaces[sharedNdx]) ? [] : vertFaces[sharedNdx]
        vertFaces[sharedNdx].push(Math.floor(i / 3));
    }
    return vertFaces;
}

function generateNormals(arrays) {
    const maxAngle = degToRad(2 * Math.PI / 180);
    const positions = arrays.positions;
    const texcoords = arrays.texcoords;
    const indices = arrays.indices;
    const faceNormals = generateFaceNormals(indices, positions);
    const vertIndices = generateVertIndices(positions);
    const vertFaces = generateVertFaces(vertIndices, indices);
    const newPositions = [];
    const newTexcoords = [];
    const newNormals = [];
    const newVertIndices = [];
    let tempVerts = {};
    let tempVertNdx = 0;
    const maxAngleCos = Math.cos(maxAngle);
    for (let i = 0; i < (indices.length); ++i) {
        const thisFaceNormal = faceNormals[Math.floor(i / 3)];
        const ndx = indices[i];
        const norm = [0, 0, 0];
        vertFaces[vertIndices[ndx]].forEach(faceNdx => {
            const otherFaceNormal = faceNormals[faceNdx];
            if (m4.dot(thisFaceNormal, otherFaceNormal) > maxAngleCos) {
                norm[0] += otherFaceNormal[0];
                norm[1] += otherFaceNormal[1];
                norm[2] += otherFaceNormal[2];
            }
        });
        m4.normalize(norm, norm);
        const poffset = ndx * 3;
        const toffset = ndx * 2;
        const x = positions[poffset + 0];
        const y = positions[poffset + 1];
        const z = positions[poffset + 2];
        const nx = norm[0];
        const ny = norm[1];
        const nz = norm[2];
        const u = texcoords[toffset + 0]
        const v = texcoords[toffset + 1]
        const vertId = `${x},${y},${z},${nx},${ny},${nz},${u},${v}`;
        if (tempVerts[vertId] === undefined) {
            tempVerts[vertId] = tempVertNdx++;
            newPositions.push(x, y, z);
            newNormals.push(nx, ny, nz);
            newTexcoords.push(u, v);
        }
        newVertIndices.push( tempVerts[vertId] );
    }
    return {
        position: newPositions,
        texcoord: newTexcoords,
        normal: newNormals,
        indices: newVertIndices,
    };
}


function generateArrays(imgData, max, min, width, height) {
    const scaleFactor =  99;
    const heightFactor = 50;
    let n = (i) => Math.floor((i - min) / (max - min) * scaleFactor) * heightFactor/scaleFactor;
    const skip_vertices = [];
    const positions = [];
    const texcoords = [];
    const indices = [];
    const cellsAcross = width - 1;
    const cellsDeep = height - 1;
    for (let z = 0; z < cellsDeep; z++) {
        for (let x = 0; x < cellsAcross; x++) {
            if (skip_vertices[`${x},${z}`] !== undefined) {
                continue;
            } 
            const ndx = (positions.length / 3);
            const base0 = z * width + x;
            const base1 = base0 + width;
            const base2 = base1 + width;
            const h00 = n(imgData[(base0 + 0)]); 
            const h01 = n(imgData[(base0 + 1)]); 
            const h10 = n(imgData[(base1 + 0)]); 
            const h11 = n(imgData[(base1 + 1)]); 
            const x0 = x;
            const z0 = z;
            let x1 = x + 1;
            let z1 = z + 1;
            const u0 = x / cellsAcross;
            const v0 = z / cellsDeep;
            let u1 = x1 / cellsAcross;
            let v1 = z1 / cellsDeep;
            if (x < cellsAcross - 2 && z < cellsDeep - 1) {
                const h02 = n(imgData[(base0 + 2)]); 
                const h12 = n(imgData[(base1 + 2)]); 
                const h20 = n(imgData[(base2 + 0)]); 
                const h21 = n(imgData[(base2 + 1)]); 
                const h22 = n(imgData[(base2 + 2)]); 
                if (h00 == h02 && h02 == h22 && h20 == h02) {
                    skip_vertices[`${x+1},${z+1}`] = 1;
                    skip_vertices[`${x+0},${z+1}`] = 1;
                    skip_vertices[`${x+1},${z+0}`] = 1;
                    x1 += 1;
                    z1 += 1;
                    u1 = x1 / cellsAcross;
                    v1 = z1 / cellsDeep;
                } else if (h00 == h01 && h02 == h01 && h12 == h11 && h10 == h11) {
                    x1 += 1;
                    u1 = x1 / cellsAcross;
                    skip_vertices[`${x+1},${z+0}`] = 1;
                } else if (h00 == h20 && h20 == h21 && h21 == h01) {
                    z1 += 1;
                    v1 = z1 / cellsDeep;
                    skip_vertices[`${x+0},${z+1}`] = 1;
                }
            }
            positions.push( x0, h00, z0,  x1, h01, z0,  x0, h10, z1,  x1, h11, z1 );
            texcoords.push( u0, v0,  u1, v0,  u0, v1,  u1, v1 );
            indices.push( ndx, ndx + 3, ndx + 1, ndx, ndx + 2, ndx + 3 );
        }
    }
    return { positions, texcoords, indices };
}

function draw3d(image) {
    let pixels;
    if (glob.efmap_altered || glob.efmap_pixels === null) {
        pixels = image.efpix();
        glob.efmap_altered = 0;
        glob.efmap_pixels = pixels;    
    } else {
        pixels = glob.efmap_pixels;
    }
    const max = Math.max(...pixels);
    const min = Math.min(...pixels);
    if (!glob.fake_canvas) {
        glob.fake_canvas = document.createElement('canvas').getContext('2d');
        glob.fake_canvas.canvas.width = glob.canvas_width;
        glob.fake_canvas.canvas.height = glob.canvas_height;
        for (let x = 0; x < glob.map_width; x++) {
            for (let y = 0; y < glob.map_height; y++) {
                let value = pixels[((y * glob.map_width) + x)];
                let newvalue = Math.floor((value - min) / (max - min) * 255);
                glob.fake_canvas.fillStyle = turbo[newvalue];
                glob.fake_canvas.fillRect(x, y, glob.cellSize, glob.cellSize);
            }
        }
        glob.refresh_canvas = 0;
    }
    const texture = glob.fake_canvas.getImageData(0, 0, glob.map_width, glob.map_height);
    const arrays = generateArrays(pixels, max, min, glob.map_width, glob.map_height);
    const arrays_with_normals = generateNormals(arrays);
    const gl = glob.gl;
    gl.clearColor(1, 1, 1, 1);
    gl.clear(gl.COLOR_BUFFER_BIT);
    gl.bindTexture(gl.TEXTURE_2D, gl.createTexture());
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, texture);
    gl.generateMipmap(gl.TEXTURE_2D);
    glob.delta_ef_map_y = 0;
    glob.delta_ef_map_x = 0;
    glob.three_d_zoom_level = 0;
    glob.programInfo = createProgramInfo(gl, [VERTEX_SHADER, FRAGMENT_SHADER]);
    glob.bufferInfo = createBufferInfoFromArrays(gl, arrays_with_normals)
    requestAnimationFrame(render);
}


function render() {
    try {
        const gl = glob.gl;
        glob.delta_ef_map_x = Math.min(89, glob.delta_ef_map_x);
        glob.delta_ef_map_x = Math.max(0, glob.delta_ef_map_x)
        glob.gl.canvas.width  = glob.gl.canvas.clientWidth  | 0;
        glob.gl.canvas.height = glob.gl.canvas.clientHeight | 0;
        const modelYRotationRadians = degToRad(-glob.delta_ef_map_y);
        const modelXRotationRadians = degToRad(-glob.delta_ef_map_x);
        const fieldOfViewRadians = degToRad(11.5);
        const aspect = gl.canvas.clientWidth / gl.canvas.clientHeight;
        const near = 10;
        const far = 2000;
        const projection = m4.perspective(fieldOfViewRadians, aspect, near, far);
        const cameraPosition = [0, 1010 + glob.three_d_zoom_level, 1];
        const target = [0, 0, 0];
        const up = [0, 1, 0];
        const camera = m4.lookAt(cameraPosition, target, up);
        let modelView = m4.inverse(camera);
        modelView = m4.xRotate(modelView, modelXRotationRadians);
        modelView = m4.yRotate(modelView, modelYRotationRadians);
        modelView = m4.translate(modelView, -glob.map_width/2, 0, -glob.map_height/2);
        let mv = modelView;
        gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
        gl.enable(gl.DEPTH_TEST);
        gl.enable(gl.CULL_FACE);
        gl.useProgram(glob.programInfo.program);
        setBuffersAndAttributes(gl, glob.programInfo, glob.bufferInfo);
        setUniforms(glob.programInfo, { projection, modelView: mv });
        const OFFSET = 0;
        gl.drawElements(gl.TRIANGLES, glob.bufferInfo.numElements, gl.UNSIGNED_INT, OFFSET);
        requestAnimationFrame(render);
    } catch (error) {
        if (glob.three_dimension_view) {
            glob.three_dimension_view = 0;
        }
    }
}

export { draw3d };

