(window.webpackJsonp=window.webpackJsonp||[]).push([[15],{17:function(t,e,r){"use strict";(function(t){r.d(e,"a",(function(){return _})),r.d(e,"b",(function(){return y}));var n=r(18);let i=new("undefined"==typeof TextDecoder?(0,t.require)("util").TextDecoder:TextDecoder)("utf-8",{ignoreBOM:!0,fatal:!0});i.decode();let s=null;function u(){return null!==s&&s.buffer===n.C.buffer||(s=new Uint8Array(n.C.buffer)),s}let o=null;function c(){return null!==o&&o.buffer===n.C.buffer||(o=new Float64Array(n.C.buffer)),o}let l=0;function a(t,e){const r=e(8*t.length);return c().set(t,r/8),l=t.length,r}function f(t,e){const r=e(1*t.length);return u().set(t,r/1),l=t.length,r}let p=null;function d(){return null!==p&&p.buffer===n.C.buffer||(p=new Int32Array(n.C.buffer)),p}function h(t,e){return c().subarray(t/8,t/8+e)}function b(t,e){return u().subarray(t/1,t/1+e)}class _{static __wrap(t){const e=Object.create(_.prototype);return e.ptr=t,e}__destroy_into_raw(){const t=this.ptr;return this.ptr=0,t}free(){const t=this.__destroy_into_raw();n.a(t)}constructor(t,e,r,i,s,u,o,c){const p=a(s,n.d),d=l,h=f(u,n.d),b=l,y=n.n(t,e,r,i,p,d,h,b,o,c);return _.__wrap(y)}efpix(){try{const i=n.b(-16);n.g(i,this.ptr);var t=d()[i/4+0],e=d()[i/4+1],r=h(t,e).slice();return n.c(t,8*e),r}finally{n.b(16)}}electrode_pixels(){try{const i=n.b(-16);n.h(i,this.ptr);var t=d()[i/4+0],e=d()[i/4+1],r=b(t,e).slice();return n.c(t,1*e),r}finally{n.b(16)}}timeline_length(){return n.r(this.ptr)>>>0}sum(){return n.q(this.ptr)>>>0}width(){return n.B(this.ptr)>>>0}height(){return n.m(this.ptr)>>>0}scale(){return n.p(this.ptr)>>>0}update_voltages(t){const e=a(t,n.d),r=l;n.A(this.ptr,e,r)}update_solids(t){const e=f(t,n.d),r=l;n.y(this.ptr,e,r)}update_offsets(t){const e=a(t,n.d),r=l;n.w(this.ptr,e,r)}update_magnets(t){const e=a(t,n.d),r=l;n.u(this.ptr,e,r)}update_pressure(t){n.x(this.ptr,t)}update_gas_mass(t){n.t(this.ptr,t)}update_fdm_threshold(t){n.s(this.ptr,t)}update_max_time(t){n.v(this.ptr,t)}update_timeline(t,e,r,i){const s=a(t,n.d),u=l,o=a(e,n.d),c=l,f=a(r,n.d),p=l,d=a(i,n.d),h=l;n.z(this.ptr,s,u,o,c,f,p,d,h)}clear(){n.f(this.ptr)}getef(t,e){return n.l(this.ptr,t,e)}save_simion_pa(){try{const i=n.b(-16);n.o(i,this.ptr);var t=d()[i/4+0],e=d()[i/4+1],r=b(t,e).slice();return n.c(t,1*e),r}finally{n.b(16)}}brush(t,e,r){n.e(this.ptr,t,e,r)}generate_electrodes(){n.k(this.ptr)}generate_electric_fields(){n.j(this.ptr)}fly_ion(t){try{const s=n.b(-16),u=a(t,n.d),o=l;n.i(s,this.ptr,u,o);var e=d()[s/4+0],r=d()[s/4+1],i=h(e,r).slice();return n.c(e,8*r),i}finally{n.b(16)}}}function y(t,e){throw new Error((r=t,n=e,i.decode(u().subarray(r,r+n))));var r,n}}).call(this,r(19)(t))},18:function(t,e,r){"use strict";var n=r.w[t.i];t.exports=n;r(17);n.D()},19:function(t,e){t.exports=function(t){if(!t.webpackPolyfill){var e=Object.create(t);e.children||(e.children=[]),Object.defineProperty(e,"loaded",{enumerable:!0,get:function(){return e.l}}),Object.defineProperty(e,"id",{enumerable:!0,get:function(){return e.i}}),Object.defineProperty(e,"exports",{enumerable:!0}),e.webpackPolyfill=1}return e}},20:function(t,e,r){"use strict";r.r(e);var n=r(17);r.d(e,"Environment",(function(){return n.a})),r.d(e,"__wbindgen_throw",(function(){return n.b}))}}]);