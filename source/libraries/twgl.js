
/*
 * Copyright 2019 Gregg Tavares
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

const defaultShaderType = [
    "VERTEX_SHADER",
    "FRAGMENT_SHADER",
];

const spaceRE = /^[ \t]*\n/;
const tokenRE = /(\.|\[|]|\w+)/g;
const texcoordRE = /coord|texture/i;

const STATIC_DRAW                    = 0x88e4;
const ELEMENT_ARRAY_BUFFER           = 0x8893;
const COMPILE_STATUS                 = 0x8b81;
const LINK_STATUS                    = 0x8b82;
const ACTIVE_UNIFORM_BLOCKS          = 0x8a36;
const ACTIVE_UNIFORMS                = 0x8b86;
const FLOAT_MAT4                     = 0x8B5C;
const ACTIVE_ATTRIBUTES              = 0x8b89;
const FLOAT_VEC4                     = 0x8B52;
const FLOAT_VEC3                     = 0x8B51;
const FLOAT_VEC2                     = 0x8B50;
const TRANSFORM_FEEDBACK_VARYINGS    = 0x8c83;
const ARRAY_BUFFER                   = 0x8892;
const BYTE                           = 0x1400;
const UNSIGNED_BYTE                  = 0x1401;
const SHORT                          = 0x1402;
const UNSIGNED_SHORT                 = 0x1403;
const INT                            = 0x1404;
const UNSIGNED_INT                   = 0x1405;
const FLOAT                          = 0x1406;
const SAMPLER_2D                     = 0x8B5E;
const TEXTURE_2D                     = 0x0DE1;


const typeMap = {};

const defaults = {
    attribPrefix: "",
};

function floatMat4Setter(gl, location) {
    return function(v) {
        gl.uniformMatrix4fv(location, false, v);
    };
}

function floatVec4Setter(gl, location) {
    return function(v) {
        gl.uniform4fv(location, v);
    };
}

function floatVec3Setter(gl, location) {
    return function(v) {
        gl.uniform3fv(location, v);
    };
}

function floatVec2Setter(gl, location) {
    return function(v) {
        gl.uniform2fv(location, v);
    };
}

function matAttribSetter(gl, index, typeInfo) {
    const defaultSize = typeInfo.size;
    const count = typeInfo.count;
  
    return function(b) {
      gl.bindBuffer(ARRAY_BUFFER, b.buffer);
      const numComponents = b.size || b.numComponents || defaultSize;
      const size = numComponents / count;
      const type = b.type || FLOAT;
      const typeInfo = typeMap[type];
      const stride = typeInfo.size * numComponents;
      const normalize = b.normalize || false;
      const offset = b.offset || 0;
      const rowOffset = stride / count;
      for (let i = 0; i < count; ++i) {
        gl.enableVertexAttribArray(index + i);
        gl.vertexAttribPointer(
            index + i, size, type, normalize, stride, offset + rowOffset * i);
        if (b.divisor !== undefined) {
          gl.vertexAttribDivisor(index + i, b.divisor);
        }
      }
    };
}

function floatAttribSetter(gl, index) {
    return function(b) {
      if (b.value) {
        gl.disableVertexAttribArray(index);
        switch (b.value.length) {
          case 4:
            gl.vertexAttrib4fv(index, b.value);
            break;
          case 3:
            gl.vertexAttrib3fv(index, b.value);
            break;
          case 2:
            gl.vertexAttrib2fv(index, b.value);
            break;
          case 1:
            gl.vertexAttrib1fv(index, b.value);
            break;
          default:
            throw new Error('the length of a float constant value must be between 1 and 4!');
        }
      } else {
        gl.bindBuffer(ARRAY_BUFFER, b.buffer);
        gl.enableVertexAttribArray(index);
        gl.vertexAttribPointer(
            index, b.numComponents || b.size, b.type || FLOAT, b.normalize || false, b.stride || 0, b.offset || 0);
        if (b.divisor !== undefined) {
          gl.vertexAttribDivisor(index, b.divisor);
        }
      }
    };
}

function samplerSetter(gl, type, unit, location) {
  var bindPoint = typeMap[type].bindPoint;
  return function(texture) {
    gl.uniform1i(location, unit);
    gl.activeTexture(gl.TEXTURE0 + unit);
    gl.bindTexture(bindPoint, texture);
  };
}

function samplerArraySetter(gl, type, unit, location, size) {
    var bindPoint = typeMap[type].bindPoint;
    var units = new Int32Array(size);
    for (var ii = 0; ii < size; ++ii) {
        units[ii] = unit + ii;
    }

    return function(textures) {
        gl.uniform1iv(location, units);
        textures.forEach(function(texture, index) {
        gl.activeTexture(gl.TEXTURE0 + units[index]);
        gl.bindTexture(bindPoint, texture);
        });
    };
}

const attrTypeMap = {};
attrTypeMap[FLOAT_VEC2]        = { size:  8, setter: floatAttribSetter, };
attrTypeMap[FLOAT_VEC3]        = { size: 12, setter: floatAttribSetter, };
attrTypeMap[FLOAT_VEC4]        = { size: 16, setter: floatAttribSetter, };
attrTypeMap[FLOAT_MAT4]        = { size: 16, setter: matAttribSetter,   count: 4, };

typeMap[FLOAT_VEC2]                    = { Type: Float32Array, size:  8, setter: floatVec2Setter,  cols: 2, };
typeMap[FLOAT_VEC3]                    = { Type: Float32Array, size: 12, setter: floatVec3Setter,  cols: 3, };
typeMap[FLOAT_VEC4]                    = { Type: Float32Array, size: 16, setter: floatVec4Setter,  cols: 4, };
typeMap[FLOAT_MAT4]                    = { Type: Float32Array, size: 64, setter: floatMat4Setter,  rows: 4, cols: 4, };
typeMap[SAMPLER_2D]                    = { Type: null,         size:  0, setter: samplerSetter,    arraySetter: samplerArraySetter, bindPoint: TEXTURE_2D,       };

function createAttributeSetters(gl, program) {
    const attribSetters = {};
    const numAttribs = gl.getProgramParameter(program, ACTIVE_ATTRIBUTES);
    for (let ii = 0; ii < numAttribs; ++ii) {
      const attribInfo = gl.getActiveAttrib(program, ii);
      if (isBuiltIn(attribInfo)) {
        continue;
      }
      const index = gl.getAttribLocation(program, attribInfo.name);
      const typeInfo = attrTypeMap[attribInfo.type];
      const setter = typeInfo.setter(gl, index, typeInfo);
      setter.location = index;
      attribSetters[attribInfo.name] = setter;
    }
    return attribSetters;
}
  

function addSetterToUniformTree(fullPath, setter, node, uniformSetters) {
    const tokens = fullPath.split(tokenRE).filter(s => s !== '');
    let tokenNdx = 0;
    let path = '';
    for (;;) {
        const token = tokens[tokenNdx++];  // has to be name or number
        path += token;
        const isArrayIndex = token[0] >= '0' && token[0] <= '9'; // isDigit
        const accessor = isArrayIndex ? parseInt(token) : token;
        if (isArrayIndex) {
            path += tokens[tokenNdx++];  // skip ']'
        }
        const isLastToken = tokenNdx === tokens.length;
        if (isLastToken) {
            node[accessor] = setter;
            break;
        } else {
            const token = tokens[tokenNdx++];  // has to be . or [
            const isArray = token === '[';
            const child = node[accessor] || (isArray ? [] : {});
            node[accessor] = child;
            node = child;
            uniformSetters[path] = uniformSetters[path] || function(node) {
            return function(value) {
                setUniformTree(node, value);
            };
            }(child);
            path += token;
        }
    }
}

function error(...args) { console.error(...args); }



function isBuiltIn(info) {
    const name = info.name;
    return name.startsWith("gl_") || name.startsWith("webgl_");
}

function createTransformFeedbackInfo(gl, program) {
    const info = {};
    const numVaryings = gl.getProgramParameter(program, TRANSFORM_FEEDBACK_VARYINGS);
    for (let ii = 0; ii < numVaryings; ++ii) {
      const varying = gl.getTransformFeedbackVarying(program, ii);
      info[varying.name] = {
        index: ii,
        type: varying.type,
        size: varying.size,
      };
    }
    return info;
}

function getNormalizationForTypedArray(typedArray) {
    if (typedArray instanceof Int8Array)    { return true; }  // eslint-disable-line
    if (typedArray instanceof Uint8Array)   { return true; }  // eslint-disable-line
    return false;
}

function createUniformSetters(gl, program) {
    let textureUnit = 0;
    function createUniformSetter(program, uniformInfo, location) {
      const isArray = uniformInfo.name.endsWith("[0]");
      const type = uniformInfo.type;
      const typeInfo = typeMap[type];
      if (!typeInfo) {
        throw new Error(`unknown type: 0x${type.toString(16)}`); // we should never get here.
      }
      let setter;
      if (typeInfo.bindPoint) {
        const unit = textureUnit;
        textureUnit += uniformInfo.size;
        if (isArray) {
          setter = typeInfo.arraySetter(gl, type, unit, location, uniformInfo.size);
        } else {
          setter = typeInfo.setter(gl, type, unit, location, uniformInfo.size);
        }
      } else {
        if (typeInfo.arraySetter && isArray) {
          setter = typeInfo.arraySetter(gl, location);
        } else {
          setter = typeInfo.setter(gl, location);
        }
      }
      setter.location = location;
      return setter;
    }
  
    const uniformSetters = {};
    const uniformTree = {};
    const numUniforms = gl.getProgramParameter(program, ACTIVE_UNIFORMS);
  
    for (let ii = 0; ii < numUniforms; ++ii) {
      const uniformInfo = gl.getActiveUniform(program, ii);
      if (isBuiltIn(uniformInfo)) {
        continue;
      }
      let name = uniformInfo.name;
      if (name.endsWith("[0]")) {
        name = name.substr(0, name.length - 3);
      }
      const location = gl.getUniformLocation(program, uniformInfo.name);
      if (location) {
        const setter = createUniformSetter(program, uniformInfo, location);
        uniformSetters[name] = setter;
        addSetterToUniformTree(name, setter, uniformTree, uniformSetters);
      }
    }
    return uniformSetters;
}

function createProgramInfoFromProgram(gl, program) {
    const uniformSetters = createUniformSetters(gl, program);
    const attribSetters = createAttributeSetters(gl, program);
    const programInfo = {
      program,
      uniformSetters,
      attribSetters,
    };
    if (!!gl.texStorage2D) {
      programInfo.uniformBlockSpec = createUniformBlockSpecFromProgram(gl, program);
      programInfo.transformFeedbackInfo = createTransformFeedbackInfo(gl, program);
    }
    return programInfo;
}

function createProgramInfo(gl, shaderSources, opt_attribs, opt_locations, opt_errorCallback) {
  const progOptions = getProgramOptions(opt_attribs, opt_locations, opt_errorCallback);
  let good = true;
  shaderSources = shaderSources.map(function(source) {
    if (source.indexOf("\n") < 0) {
      const script = getElementById(source);
      if (!script) {
        progOptions.errorCallback("no element with id: " + source);
        good = false;
      } else {
        source = script.text;
      }
    }
    return source;
  });
  if (!good) {
    return null;
  }
  const program = createProgramFromSources(gl, shaderSources, progOptions);
  if (!program) {
    return null;
  }
  return createProgramInfoFromProgram(gl, program);
}

function loadShader(gl, shaderSource, shaderType, opt_errorCallback) {
    const errFn = opt_errorCallback || error;
    const shader = gl.createShader(shaderType);
    let lineOffset = 0;
    if (spaceRE.test(shaderSource)) {
      lineOffset = 1;
      shaderSource = shaderSource.replace(spaceRE, '');
    }
    gl.shaderSource(shader, shaderSource);
    gl.compileShader(shader);
    const compiled = gl.getShaderParameter(shader, COMPILE_STATUS);
    if (!compiled) {
      const lastError = gl.getShaderInfoLog(shader);
      errFn(`${addLineNumbersWithError(shaderSource, lastError, lineOffset)}\nError compiling ${glEnumToString(gl, shaderType)}: ${lastError}`);
      gl.deleteShader(shader);
      return null;
    }
    return shader;
}

function createProgramFromSources(gl, shaderSources, opt_attribs, opt_locations, opt_errorCallback) {
  const progOptions = getProgramOptions(opt_attribs, opt_locations, opt_errorCallback);
  const shaders = [];
  for (let ii = 0; ii < shaderSources.length; ++ii) {
    const shader = loadShader(
        gl, shaderSources[ii], gl[defaultShaderType[ii]], progOptions.errorCallback);
    if (!shader) {
      return null;
    }
    shaders.push(shader);
  }
  return createProgram(gl, shaders, progOptions);
}

function getProgramOptions(opt_attribs, opt_locations, opt_errorCallback) {
    let transformFeedbackVaryings;
    let transformFeedbackMode;
    if (typeof opt_locations === 'function') {
      opt_errorCallback = opt_locations;
      opt_locations = undefined;
    }
    if (typeof opt_attribs === 'function') {
      opt_errorCallback = opt_attribs;
      opt_attribs = undefined;
    } else if (opt_attribs && !Array.isArray(opt_attribs)) {
      // If we have an errorCallback we can just return this object
      // Otherwise we need to construct one with default errorCallback
      if (opt_attribs.errorCallback) {
        return opt_attribs;
      }
      const opt = opt_attribs;
      opt_errorCallback = opt.errorCallback;
      opt_attribs = opt.attribLocations;
      transformFeedbackVaryings = opt.transformFeedbackVaryings;
      transformFeedbackMode = opt.transformFeedbackMode;
    }
  
    const options = {
      errorCallback: opt_errorCallback || error,
      transformFeedbackVaryings: transformFeedbackVaryings,
      transformFeedbackMode: transformFeedbackMode,
    };
  
    if (opt_attribs) {
      let attribLocations = {};
      if (Array.isArray(opt_attribs)) {
        opt_attribs.forEach(function(attrib,  ndx) {
          attribLocations[attrib] = opt_locations ? opt_locations[ndx] : ndx;
        });
      } else {
        attribLocations = opt_attribs;
      }
      options.attribLocations = attribLocations;
    }
  
    return options;
}

function createProgram(
    gl, shaders, opt_attribs, opt_locations, opt_errorCallback) {
    const progOptions = getProgramOptions(opt_attribs, opt_locations, opt_errorCallback);
    const realShaders = [];
    const newShaders = [];
    for (let ndx = 0; ndx < shaders.length; ++ndx) {
        let shader = shaders[ndx];
        if (typeof (shader) === 'string') {
        const elem = getElementById(shader);
        const src = elem ? elem.text : shader;
        let type = gl[defaultShaderType[ndx]];
        if (elem && elem.type) {
            type = getShaderTypeFromScriptType(gl, elem.type) || type;
        }
        shader = loadShader(gl, src, type, progOptions.errorCallback);
        newShaders.push(shader);
        }
        if (typeof WebGLShader !== 'undefined' && shader instanceof WebGLShader) {
            realShaders.push(shader);
        }
    }

    if (realShaders.length !== shaders.length) {
        progOptions.errorCallback("not enough shaders for program");
        deleteShaders(gl, newShaders);
        return null;
    }

    const program = gl.createProgram();
    realShaders.forEach(function(shader) {
        gl.attachShader(program, shader);
    });
    if (progOptions.attribLocations) {
        Object.keys(progOptions.attribLocations).forEach(function(attrib) {
        gl.bindAttribLocation(program, progOptions.attribLocations[attrib], attrib);
        });
    }
    let varyings = progOptions.transformFeedbackVaryings;
    if (varyings) {
        if (varyings.attribs) {
            varyings = varyings.attribs;
        }
        if (!Array.isArray(varyings)) {
            varyings = Object.keys(varyings);
        }
        gl.transformFeedbackVaryings(program, varyings, progOptions.transformFeedbackMode || SEPARATE_ATTRIBS);
    }
    gl.linkProgram(program);

    // Check the link status
    const linked = gl.getProgramParameter(program, LINK_STATUS);
    if (!linked) {
        // something went wrong with the link
        const lastError = gl.getProgramInfoLog(program);
        progOptions.errorCallback(`${
        realShaders.map(shader => {
            const src = addLineNumbersWithError(gl.getShaderSource(shader), '', 0);
            const type = gl.getShaderParameter(shader, gl.SHADER_TYPE);
            return `${glEnumToString(gl, type)}\n${src}}`;
        }).join('\n')
        }\nError in program linking: ${lastError}`);

        gl.deleteProgram(program);
        deleteShaders(gl, newShaders);
        return null;
    }
    return program;
}

function createBufferInfoFromArrays(gl, arrays, srcBufferInfo) {
    const newAttribs = createAttribsFromArrays(gl, arrays);
    const bufferInfo = Object.assign({}, srcBufferInfo ? srcBufferInfo : {});
    bufferInfo.attribs = Object.assign({}, srcBufferInfo ? srcBufferInfo.attribs : {}, newAttribs);
    const indices = arrays.indices;
    if (indices) {
      const newIndices = makeTypedArray(indices, "indices");
      bufferInfo.indices = createBufferFromTypedArray(gl, newIndices, ELEMENT_ARRAY_BUFFER);
      bufferInfo.numElements = newIndices.length;
      bufferInfo.elementType = getGLTypeForTypedArray(newIndices);
    } else if (!bufferInfo.numElements) {
      bufferInfo.numElements = getNumElementsFromAttributes(gl, bufferInfo.attribs);
    }
    return bufferInfo;
}

function createBufferFromTypedArray(gl, typedArray, type, drawType) {
    type = type || ARRAY_BUFFER;
    const buffer = gl.createBuffer();
    setBufferFromTypedArray(gl, type, buffer, typedArray, drawType);
    return buffer;
}

function setBufferFromTypedArray(gl, type, buffer, array, drawType) {
    gl.bindBuffer(type, buffer);
    gl.bufferData(type, array, drawType || STATIC_DRAW);
}

function setBuffersAndAttributes(gl, programInfo, buffers) {
    if (buffers.vertexArrayObject) {
        gl.bindVertexArray(buffers.vertexArrayObject);
    } else {
        setAttributes(programInfo.attribSetters || programInfo, buffers.attribs);
        if (buffers.indices) {
            gl.bindBuffer(ELEMENT_ARRAY_BUFFER, buffers.indices);
        }
    }
}

function setAttributes(setters, buffers) {
    for (const name in buffers) {
        const setter = setters[name];
        if (setter) {
            setter(buffers[name]);
        }
    }
}

function setUniforms(setters, ...args) {  // eslint-disable-line
    const actualSetters = setters.uniformSetters || setters;
    const numArgs = args.length;
    for (let aNdx = 0; aNdx < numArgs; ++aNdx) {
        const values = args[aNdx];
        if (Array.isArray(values)) {
            const numValues = values.length;
            for (let ii = 0; ii < numValues; ++ii) {
                setUniforms(actualSetters, values[ii]);
            }
        } else {
            for (const name in values) {
                const setter = actualSetters[name];
                if (setter) {
                    setter(values[name]);
                }
            }
        }
    }
}


function createAttribsFromArrays(gl, arrays) {
    const attribs = {};
    Object.keys(arrays).forEach(function(arrayName) {
      if (!(arrayName === "indices")) {
        const array = arrays[arrayName];
        const attribName = array.attrib || array.name || array.attribName || (defaults.attribPrefix + arrayName);
        if (array.value) {
          if (!Array.isArray(array.value) && !isArrayBuffer(array.value)) {
            throw new Error('array.value is not array or typedarray');
          }
          attribs[attribName] = {
            value: array.value,
          };
        } else {
          let buffer;
          let type;
          let normalization;
          let numComponents;
          if (array.buffer && array.buffer instanceof WebGLBuffer) {
            buffer = array.buffer;
            numComponents = array.numComponents || array.size;
            type = array.type;
            normalization = array.normalize;
          } else if (typeof array === "number" || typeof array.data === "number") {
            const numValues = array.data || array;
            const arrayType = array.type || Float32Array;
            const numBytes = numValues * arrayType.BYTES_PER_ELEMENT;
            type = getGLTypeForTypedArrayType(arrayType);
            normalization = array.normalize !== undefined ? array.normalize : getNormalizationForTypedArrayType(arrayType);
            numComponents = array.numComponents || array.size || guessNumComponentsFromName(arrayName, numValues);
            buffer = gl.createBuffer();
            gl.bindBuffer(ARRAY_BUFFER, buffer);
            gl.bufferData(ARRAY_BUFFER, numBytes, array.drawType || STATIC_DRAW);
          } else {
            const typedArray = makeTypedArray(array, arrayName);
            buffer = createBufferFromTypedArray(gl, typedArray, undefined, array.drawType);
            type = getGLTypeForTypedArray(typedArray);
            normalization = array.normalize !== undefined ? array.normalize : getNormalizationForTypedArray(typedArray);
            numComponents = getNumComponents(array, arrayName);
          }
          attribs[attribName] = {
            buffer:        buffer,
            numComponents: numComponents,
            type:          type,
            normalize:     normalization,
            stride:        array.stride || 0,
            offset:        array.offset || 0,
            divisor:       array.divisor === undefined ? undefined : array.divisor,
            drawType:      array.drawType,
          };
        }
      }
    });
    gl.bindBuffer(ARRAY_BUFFER, null);
    return attribs;
}

function getGLTypeForTypedArray(typedArray) {
    if (typedArray instanceof Int8Array)         { return BYTE; }           // eslint-disable-line
    if (typedArray instanceof Uint8Array)        { return UNSIGNED_BYTE; }  // eslint-disable-line
    if (typedArray instanceof Uint8ClampedArray) { return UNSIGNED_BYTE; }  // eslint-disable-line
    if (typedArray instanceof Int16Array)        { return SHORT; }          // eslint-disable-line
    if (typedArray instanceof Uint16Array)       { return UNSIGNED_SHORT; } // eslint-disable-line
    if (typedArray instanceof Int32Array)        { return INT; }            // eslint-disable-line
    if (typedArray instanceof Uint32Array)       { return UNSIGNED_INT; }   // eslint-disable-line
    if (typedArray instanceof Float32Array)      { return FLOAT; }          // eslint-disable-line
    throw new Error('unsupported typed array type');
}

const isArrayBuffer = typeof SharedArrayBuffer !== 'undefined'
  ? function isArrayBufferOrSharedArrayBuffer(a) {
    return a && a.buffer && (a.buffer instanceof ArrayBuffer || a.buffer instanceof SharedArrayBuffer);
  }
  : function isArrayBuffer(a) {
    return a && a.buffer && a.buffer instanceof ArrayBuffer;
  };

function makeTypedArray(array, name) {
    if (isArrayBuffer(array)) {
      return array;
    }
  
    if (isArrayBuffer(array.data)) {
      return array.data;
    }
  
    if (Array.isArray(array)) {
      array = {
        data: array,
      };
    }
  
    let Type = array.type;
    if (!Type) {
      if (name === "indices") {
        Type = Uint32Array;
      } else {
        Type = Float32Array;
      }
    }
    return new Type(array.data);
}

function createUniformBlockSpecFromProgram(gl, program) {
    const numUniforms = gl.getProgramParameter(program, ACTIVE_UNIFORMS);
    const uniformData = [];
    const uniformIndices = [];
  
    for (let ii = 0; ii < numUniforms; ++ii) {
      uniformIndices.push(ii);
      uniformData.push({});
      const uniformInfo = gl.getActiveUniform(program, ii);
      uniformData[ii].name = uniformInfo.name;
    }
  
    [
      [ "UNIFORM_TYPE", "type" ],
      [ "UNIFORM_SIZE", "size" ],  // num elements
      [ "UNIFORM_BLOCK_INDEX", "blockNdx" ],
      [ "UNIFORM_OFFSET", "offset", ],
    ].forEach(function(pair) {
      const pname = pair[0];
      const key = pair[1];
      gl.getActiveUniforms(program, uniformIndices, gl[pname]).forEach(function(value, ndx) {
        uniformData[ndx][key] = value;
      });
    });
  
    const blockSpecs = {};
  
    const numUniformBlocks = gl.getProgramParameter(program, ACTIVE_UNIFORM_BLOCKS);
    for (let ii = 0; ii < numUniformBlocks; ++ii) {
      const name = gl.getActiveUniformBlockName(program, ii);
      const blockSpec = {
        index: gl.getUniformBlockIndex(program, name),
        usedByVertexShader: gl.getActiveUniformBlockParameter(program, ii, UNIFORM_BLOCK_REFERENCED_BY_VERTEX_SHADER),
        usedByFragmentShader: gl.getActiveUniformBlockParameter(program, ii, UNIFORM_BLOCK_REFERENCED_BY_FRAGMENT_SHADER),
        size: gl.getActiveUniformBlockParameter(program, ii, UNIFORM_BLOCK_DATA_SIZE),
        uniformIndices: gl.getActiveUniformBlockParameter(program, ii, UNIFORM_BLOCK_ACTIVE_UNIFORM_INDICES),
      };
      blockSpec.used = blockSpec.usedByVertexShader || blockSpec.usedByFragmentShader;
      blockSpecs[name] = blockSpec;
    }
  
    return {
      blockSpecs: blockSpecs,
      uniformData: uniformData,
    };
}

const glEnumToString = (function() {
    const haveEnumsForType = {};
    const enums = {};

    function addEnums(gl) {
        const type = gl.constructor.name;
        if (!haveEnumsForType[type]) {
        for (const key in gl) {
            if (typeof gl[key] === 'number') {
            const existing = enums[gl[key]];
            enums[gl[key]] = existing ? `${existing} | ${key}` : key;
            }
        }
        haveEnumsForType[type] = true;
        }
    }

    return function glEnumToString(gl, value) {
        addEnums(gl);
        return enums[value] || (typeof value === 'number' ? `0x${value.toString(16)}` : value);
    };
    }()
);

function getNumComponents(array, arrayName) {
    return array.numComponents || array.size || guessNumComponentsFromName(arrayName, (array.length ? array:array.data).length);
}

function guessNumComponentsFromName(name, length) {
    let numComponents;
    if (texcoordRE.test(name)) {
        numComponents = 2;
    } else {
        numComponents = 3;
    }
    return numComponents;
}

export { createProgramInfo, createBufferInfoFromArrays, setUniforms, setBuffersAndAttributes }