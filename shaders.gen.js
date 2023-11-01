/* File generated with Shader Minifier 1.2
 * http://www.ctrl-alt-test.fr
 */
export let aimRay_frag = `precision highp float;uniform vec3 a0;void main(){gl_FragColor=vec4(a0,.25);}`;
export let aimRay_vert = `attribute vec3 a1;uniform mat4 a2,a3;void main(){gl_Position=a2*a3*vec4(a1,1);}`;
export let blit_frag = `precision highp float;uniform sampler2D a4;uniform vec4 a5;void main(){vec2 a=gl_FragCoord.xy/a5.xy;vec4 i=texture2D(a4,vec2(a.x,1.-a.y));gl_FragColor=mix(vec4(mix(vec3(0),i.xyz,a5.w),1),i,a5.z);}`;
export let blit_vert = `attribute vec2 a1;void main(){gl_Position=vec4(a1,0,1);}`;
export let debugGeo_frag = `precision highp float;varying vec3 a6,a7;uniform vec3 a0;void main(){gl_FragColor=vec4(a0,1);}`;
export let debugLines_frag = `precision highp float;varying float a8;void main(){gl_FragColor=a8>1.5?vec4(0,1,1,1):a8>.5?vec4(1,0,0,1):vec4(1,1,0,1);}`;
export let debugLines_vert = `attribute vec3 a1;attribute float a9;uniform mat4 aa;varying float a8;void main(){a8=a9;vec4 e=aa*vec4(a1,1);e.z-=.001;gl_Position=e;}`;
export let debugRay_frag = `precision highp float;varying vec3 ab;uniform vec3 a0;void main(){float o=.1*(ab.x+ab.y+ab.z);gl_FragColor=vec4((fract(o)>.5?1.:0.)*a0,1);}`;
export let debugRay_vert = `attribute float ac;uniform mat4 aa;uniform vec3 ad,ae;varying vec3 ab;void main(){vec3 e=mix(ad,ae,ac);ab=e;gl_Position=aa*vec4(e,1);}`;
export let main_frag = `precision highp float;varying vec3 af,a6;varying float a8;uniform sampler2D a4[7];uniform vec3 ag;void main(){vec3 v=floor(fract(af*.005+.01)*64.);vec2 a=(v.xy+vec2(mod(v.z,8.),floor(v.z/8.))*64.+vec2(.5))/512.;vec4 i;int t=int(a8+.5);if(t==0)i=texture2D(a4[0],a);if(t==1)i=texture2D(a4[1],a);if(t==2)i=texture2D(a4[2],a);if(t==3)i=texture2D(a4[3],a);if(t==4)i=texture2D(a4[4],a);if(t==5)i=texture2D(a4[5],a);if(t==6)i=texture2D(a4[6],a);float r=dot(a6,normalize(vec3(5,1,5)));vec3 c=mix(mix(vec3(.9,.42,.44)+.1,vec3(1),.3),mix(vec3(1,.86,.39),vec3(1),.3),.5+.5*r);float n=.75+.25*max(0.,r),f=length(af.xz-ag.xz);if(f<20.&&af.y<ag.y)n*=1.-clamp(.5-.1*(f-5.),0.,.4);gl_FragColor=vec4(i.xyz*c*n,1);}`;
export let main_vert = `attribute vec3 a1,ah;attribute float a9;uniform mat4 a2,a3;varying vec3 af,a6;varying float a8;void main(){af=a1,a6=normalize((a3*vec4(ah,0)).xyz),a8=a9,gl_Position=a2*a3*vec4(a1,1);}`;
export let sky_frag = `precision highp float;varying vec3 af;void main(){vec3 m=normalize(af);float l=dot(normalize(vec3(5,1,5)),m),g=max(0.,500.*l-495.);vec3 x=mix(vec3(.13,.16,.31),vec3(.9,.42,.44),1.-2.*m.y),u=mix(mix(x,vec3(1,.86,.39),clamp(exp(l-1.)-.5*m.y,0.,1.)),vec3(1),g);gl_FragColor=vec4(u,1);}`;
export let sky_vert = `attribute vec3 a1,ah;attribute float a9;uniform mat4 aa;varying vec3 af;void main(){af=a1,gl_Position=(aa*vec4(a1,1)).xyww;}`;
