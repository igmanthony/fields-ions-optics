(window.webpackJsonp=window.webpackJsonp||[]).push([[12],{14:function(t,r,e){"use strict";(function(t){e.d(r,"a",(function(){return _})),e.d(r,"b",(function(){return y}));var n=e(15);let i=new("undefined"==typeof TextDecoder?(0,t.require)("util").TextDecoder:TextDecoder)("utf-8",{ignoreBOM:!0,fatal:!0});i.decode();let u=null;function s(){return null!==u&&u.buffer===n.C.buffer||(u=new Uint8Array(n.C.buffer)),u}let a=null;function c(){return null!==a&&a.buffer===n.C.buffer||(a=new Float64Array(n.C.buffer)),a}let o=0;function l(t,r){const e=r(8*t.length);return c().set(t,e/8),o=t.length,e}function f(t,r){const e=r(1*t.length);return s().set(t,e/1),o=t.length,e}let p=null;function d(){return null!==p&&p.buffer===n.C.buffer||(p=new Int32Array(n.C.buffer)),p}function h(t,r){return c().subarray(t/8,t/8+r)}function b(t,r){return s().subarray(t/1,t/1+r)}class _{static __wrap(t){const r=Object.create(_.prototype);return r.ptr=t,r}__destroy_into_raw(){const t=this.ptr;return this.ptr=0,t}free(){const t=this.__destroy_into_raw();n.a(t)}constructor(t,r,e,i,u,s,a,c){var p=l(u,n.d),d=o,h=f(s,n.d),b=o,y=n.n(t,r,e,i,p,d,h,b,a,c);return _.__wrap(y)}efpix(){try{const i=n.b(-16);n.g(i,this.ptr);var t=d()[i/4+0],r=d()[i/4+1],e=h(t,r).slice();return n.c(t,8*r),e}finally{n.b(16)}}electrode_pixels(){try{const i=n.b(-16);n.h(i,this.ptr);var t=d()[i/4+0],r=d()[i/4+1],e=b(t,r).slice();return n.c(t,1*r),e}finally{n.b(16)}}timeline_length(){return n.s(this.ptr)>>>0}sum(){return n.r(this.ptr)>>>0}width(){return n.B(this.ptr)>>>0}height(){return n.m(this.ptr)>>>0}scale(){return n.q(this.ptr)>>>0}update_voltages(t){var r=l(t,n.d),e=o;n.A(this.ptr,r,e)}update_solids(t){var r=f(t,n.d),e=o;n.y(this.ptr,r,e)}update_offsets(t){var r=l(t,n.d),e=o;n.w(this.ptr,r,e)}update_magnets(t){var r=l(t,n.d),e=o;n.v(this.ptr,r,e)}update_pressure(t){n.x(this.ptr,t)}update_gas_mass(t){n.u(this.ptr,t)}update_fdm_threshold(t){n.t(this.ptr,t)}update_timeline(t,r,e,i){var u=l(t,n.d),s=o,a=l(r,n.d),c=o,f=l(e,n.d),p=o,d=l(i,n.d),h=o;n.z(this.ptr,u,s,a,c,f,p,d,h)}clear(){n.f(this.ptr)}getef(t,r){return n.l(this.ptr,t,r)}scalars(t){try{const u=n.b(-16);n.p(u,this.ptr,t);var r=d()[u/4+0],e=d()[u/4+1],i=h(r,e).slice();return n.c(r,8*e),i}finally{n.b(16)}}save_simion_pa(){try{const i=n.b(-16);n.o(i,this.ptr);var t=d()[i/4+0],r=d()[i/4+1],e=b(t,r).slice();return n.c(t,1*r),e}finally{n.b(16)}}brush(t,r,e){n.e(this.ptr,t,r,e)}generate_electrodes(){n.k(this.ptr)}generate_electric_fields(){n.j(this.ptr)}fly_ion(t){try{const a=n.b(-16);var r=l(t,n.d),e=o;n.i(a,this.ptr,r,e);var i=d()[a/4+0],u=d()[a/4+1],s=h(i,u).slice();return n.c(i,8*u),s}finally{n.b(16)}}}const y=function(t,r){throw new Error((e=t,n=r,i.decode(s().subarray(e,e+n))));var e,n}}).call(this,e(16)(t))},15:function(t,r,e){"use strict";var n=e.w[t.i];t.exports=n;e(14);n.D()},16:function(t,r){t.exports=function(t){if(!t.webpackPolyfill){var r=Object.create(t);r.children||(r.children=[]),Object.defineProperty(r,"loaded",{enumerable:!0,get:function(){return r.l}}),Object.defineProperty(r,"id",{enumerable:!0,get:function(){return r.i}}),Object.defineProperty(r,"exports",{enumerable:!0}),r.webpackPolyfill=1}return r}},17:function(t,r,e){"use strict";e.r(r);var n=e(14);e.d(r,"Environment",(function(){return n.a})),e.d(r,"__wbindgen_throw",(function(){return n.b}))}}]);