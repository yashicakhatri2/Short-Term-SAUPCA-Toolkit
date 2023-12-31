function getFinalDelOffset = getFinalDelOffset(E,Ecc,G,H,J2,L,RE,TA,b,cs,g,h,hs,ks,l,mu,r,sSun)
%getFinalDelOffset
%    getFinalDelOffset = getFinalDelOffset(E,Ecc,G,H,J2,L,RE,TA,B,CS,g,h,HS,KS,l,MU,R,sSun)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    27-May-2022 10:31:29

t2 = conj(E);
t3 = conj(Ecc);
t4 = conj(G);
t5 = conj(H);
t6 = conj(J2);
t7 = conj(L);
t8 = conj(RE);
t9 = conj(TA);
t10 = conj(b);
t11 = conj(cs);
t12 = conj(g);
t13 = conj(h);
t14 = conj(hs);
t15 = conj(ks);
t16 = conj(l);
t17 = conj(mu);
t18 = conj(r);
t19 = conj(sSun);
t28 = G.^2;
t29 = H.^2;
t52 = 1.0./L.^2;
t20 = cos(t3);
t21 = cos(t9);
t22 = sin(t3);
t23 = sin(t9);
t24 = cos(t12);
t25 = cos(t15);
t26 = sin(t12);
t27 = sin(t15);
t30 = t2.^6;
t31 = t3.*2.0;
t32 = t3.*6.0;
t33 = t4.^2;
t34 = t4.^3;
t36 = t5.^2;
t37 = t6.^2;
t38 = t7.^2;
t39 = t7.^3;
t40 = t7.^5;
t42 = t8.^2;
t44 = t9.*2.0;
t45 = t9.*3.0;
t46 = t11+1.0;
t47 = t12.*2.0;
t48 = t16.*6.0;
t49 = t17.^2;
t51 = 1.0./t28;
t53 = 1.0./t2.^3;
t54 = 1.0./t2.^5;
t56 = 1.0./t4;
t62 = 1.0./t7;
t67 = t11-1.0;
t68 = -t13;
t69 = -t14;
t70 = -t15;
t71 = -t16;
t73 = 1.0./t17;
t77 = 1.0./t18;
t96 = t28.*t52;
t35 = t33.^2;
t41 = t38.^3;
t43 = t42.^2;
t50 = t49.^2;
t55 = 1.0./t30;
t57 = 1.0./t33;
t58 = 1.0./t34;
t61 = t56.^7;
t63 = 1.0./t38;
t64 = 1.0./t39;
t72 = -t48;
t74 = 1.0./t49;
t75 = t73.^3;
t78 = t77.^2;
t79 = t77.^3;
t81 = cos(t31);
t82 = cos(t44);
t83 = t21.^2;
t84 = sin(t31);
t85 = t22.^2;
t86 = sin(t44);
t87 = t23.^2;
t88 = t23.^3;
t89 = cos(t47);
t91 = t9+t47;
t92 = t5.*t56;
t93 = t9+t71;
t95 = t29.*t51;
t99 = t44+t47;
t100 = t45+t47;
t110 = -t96;
t125 = t12+t13+t15+t69;
t126 = t12+t14+t15+t68;
t140 = t12+t13+t69+t70;
t141 = t12+t14+t68+t70;
t59 = 1.0./t35;
t60 = t57.^3;
t65 = t63.^2;
t66 = 1.0./t41;
t76 = 1.0./t50;
t80 = t78.^2;
t90 = t58.^3;
t94 = cos(t91);
t97 = sin(t91);
t98 = t4.*t63.*2.0;
t101 = cos(t99);
t102 = cos(t100);
t103 = sin(t99);
t104 = sin(t100);
t105 = t92+1.0;
t106 = t36.*t57;
t107 = t33.*t63;
t108 = t34.*t64;
t109 = -t95;
t111 = t92-1.0;
t113 = t33.*t64.*2.0;
t115 = t33.*t64.*3.0;
t118 = t110+1.0;
t127 = cos(t125);
t128 = cos(t126);
t129 = sin(t125);
t130 = sin(t126);
t138 = t32+t72+t84;
t148 = cos(t140);
t149 = cos(t141);
t150 = sin(t140);
t151 = sin(t141);
t152 = t4.*t39.*t74.*t78;
t112 = t106.*3.0;
t114 = t107.*4.0;
t116 = t34.*t65.*3.0;
t117 = t109+1.0;
t119 = t106.*8.0;
t120 = t106.*1.5e+1;
t121 = t106.*4.2e+1;
t122 = -t108;
t123 = -t115;
t132 = t106-1.0;
t133 = t107-1.0;
t134 = t103./2.0;
t146 = sqrt(t118);
t162 = t152-1.0;
t177 = (t46.*t56.*t128)./4.0;
t178 = (t46.*t56.*t130)./4.0;
t180 = t101.*t152;
t182 = (t56.*t67.*t127)./4.0;
t183 = (t56.*t67.*t129)./4.0;
t185 = (t5.*t46.*t57.*t128)./4.0;
t186 = (t5.*t46.*t57.*t130)./4.0;
t188 = (t46.*t56.*t148)./4.0;
t189 = t67.*t105.*t127;
t190 = t46.*t111.*t128;
t191 = (t46.*t56.*t150)./4.0;
t192 = t67.*t105.*t129;
t193 = t46.*t111.*t130;
t194 = (t5.*t57.*t67.*t127)./4.0;
t195 = (t5.*t57.*t67.*t129)./4.0;
t200 = (t56.*t67.*t149)./4.0;
t201 = (t56.*t67.*t151)./4.0;
t211 = t46.*t105.*t148;
t212 = t46.*t105.*t150;
t213 = (t5.*t46.*t57.*t148)./4.0;
t214 = (t5.*t46.*t57.*t150)./4.0;
t220 = (t5.*t57.*t67.*t149)./4.0;
t221 = (t5.*t57.*t67.*t151)./4.0;
t226 = t67.*t111.*t149;
t227 = t67.*t111.*t151;
t124 = -t116;
t131 = -t119;
t135 = t112-1.0;
t137 = t114-4.0;
t139 = t132.^2;
t142 = t120-1.0;
t144 = 1.0./t133;
t145 = sqrt(t117);
t147 = t146.^3;
t155 = conj(t146);
t165 = t98+t123;
t166 = t107+t122;
t184 = -t178;
t187 = -t183;
t196 = -t186;
t199 = -t188;
t202 = -t195;
t204 = -t200;
t206 = -t189;
t207 = -t190;
t215 = t189./4.0;
t216 = t190./4.0;
t217 = t192./4.0;
t218 = t193./4.0;
t219 = -t213;
t228 = -t212;
t231 = -t220;
t235 = t211./4.0;
t236 = t212./4.0;
t244 = t226./4.0;
t245 = t227./4.0;
t136 = t135-2.0;
t143 = t135.^2;
t153 = t139.*5.0;
t154 = conj(t145);
t156 = conj(t147);
t158 = t139.*2.1e+1;
t160 = 1.0./t155;
t163 = t21.*t155;
t164 = t23.*t155;
t168 = t113+t124;
t169 = t5.*t57.*t135.*4.8e+1;
t171 = t36.*t58.*t135.*4.8e+1;
t176 = t94.*t155;
t198 = (t97.*t155)./2.0;
t203 = (t102.*t155)./3.0;
t205 = (t104.*t155)./6.0;
t222 = -t215;
t223 = -t216;
t224 = -t218;
t239 = -t235;
t240 = -t236;
t248 = -t245;
t276 = (t102.*t152.*t155)./2.0;
t157 = t143.*4.0;
t159 = 1.0./t154;
t161 = 1.0./t156;
t167 = t163+2.0;
t172 = t19.*t24.*t25.*t154;
t173 = t19.*t24.*t27.*t154;
t174 = t19.*t25.*t26.*t154;
t175 = t19.*t26.*t27.*t154;
t179 = -t171;
t197 = t93+t164;
t208 = t131+t153;
t209 = t4.*t21.*t63.*t160;
t210 = t4.*t23.*t63.*t160;
t229 = t21.*t33.*t64.*t160;
t230 = t23.*t33.*t64.*t160;
t232 = t121+t158-2.2e+1;
t246 = t4.*t63.*t94.*t160;
t247 = t152.*t163;
t249 = t21.*t34.*t74.*t78.*t160;
t253 = (t4.*t63.*t97.*t160)./2.0;
t257 = (t4.*t63.*t102.*t160)./3.0;
t258 = (t4.*t63.*t104.*t160)./6.0;
t259 = (t33.*t64.*t97.*t160)./2.0;
t262 = t107.*(t119-t153).*-3.0;
t263 = (t33.*t64.*t104.*t160)./6.0;
t267 = t23.*t143.*t166.*8.0;
t272 = (t152.*t176)./2.0;
t273 = t86.*t143.*t166.*2.0;
t274 = (t34.*t74.*t78.*t94.*t160)./2.0;
t277 = (t34.*t74.*t78.*t102.*t160)./2.0;
t288 = t101+t176+t203;
t289 = t89.*t93.*t107.*t132.*t142.*t156.*6.0;
t294 = t134+t198+t205;
t306 = t215+t223+t239+t244;
t308 = t217+t224+t240+t245;
t170 = t167.^2;
t181 = t175.*4.0;
t225 = -t209;
t233 = t56.*t87.*t167;
t234 = t62.*t87.*t167;
t237 = t21.*t23.*t56.*t167;
t238 = t21.*t23.*t62.*t167;
t241 = -t229;
t242 = t5.*t19.*t24.*t27.*t57.*t159;
t243 = t5.*t19.*t26.*t27.*t57.*t159;
t250 = t19.*t24.*t27.*t36.*t58.*t159;
t251 = -t246;
t252 = t19.*t26.*t27.*t36.*t58.*t159;
t254 = t36.*t58.*t197.*6.0;
t255 = t23.*t56.*t97.*t167;
t260 = -t257;
t261 = t23.*t56.*t104.*t167;
t264 = (t23.*t56.*t94.*t167)./2.0;
t265 = (t23.*t62.*t94.*t167)./2.0;
t266 = t135.*t197;
t268 = (t23.*t56.*t102.*t167)./2.0;
t269 = (t23.*t62.*t102.*t167)./2.0;
t270 = t4.*t38.*t74.*t78.*t87.*t167;
t271 = t114.*t232;
t278 = t23.*t56.*t160.*t167;
t279 = t23.*t62.*t160.*t167;
t280 = (t4.*t23.*t38.*t74.*t78.*t97.*t167)./2.0;
t281 = t4.*t23.*t38.*t74.*t78.*t104.*t167.*(3.0./2.0);
t282 = t162+t247;
t283 = t157+t262;
t293 = t4.*t23.*t38.*t74.*t78.*t103.*t160.*t167.*2.0;
t296 = t36.*t58.*t294.*6.0;
t297 = t136.*t294;
t298 = t6.*t42.*t49.*t53.*t59.*t136.*t288.*(9.0./2.0);
t300 = t180+t272+t276;
t314 = t174+t217+t218+t236+t245;
t317 = t172+t215+t216+t235+t244;
t319 = t173+t217+t218+t240+t248;
t320 = t175+t222+t223+t235+t244;
t256 = -t254;
t275 = -t266;
t284 = t101.*t278;
t285 = t101.*t279;
t286 = t103.*t278.*2.0;
t287 = t225+t233;
t290 = t234+t241;
t291 = t135.*t282;
t295 = t133.*t283;
t299 = t210+t237+t278;
t301 = t230+t238+t279;
t303 = t136.*t300;
t311 = t177+t182+t199+t204+t243;
t312 = t184+t187+t191+t201+t242;
t313 = t181+t206+t207+t211+t226;
t315 = t185+t194+t219+t231+t252;
t316 = t196+t202+t214+t221+t250;
t318 = t6.*t42.*t49.*t53.*t59.*(t266-t297).*(-9.0./2.0);
t321 = t4.*t62.*t81.*t319;
t323 = t138.*t320;
t292 = -t291;
t302 = t135.*t299;
t304 = t135.*t301;
t307 = t275+t297;
t309 = t157+t271+t295;
t322 = -t321;
t327 = t251+t255+t260+t261+t286;
t328 = t253+t258+t264+t268+t284;
t329 = t259+t263+t265+t269+t285;
t305 = -t302;
t310 = t93.*t309;
t325 = t292+t303;
t330 = t136.*t328;
t331 = t136.*t329;
t334 = t322+t323;
t324 = t273+t310;
t332 = -t331;
t336 = t256+t296+t305+t330;
t337 = t6.*t42.*t49.*t53.*t58.*(t254-t296+t302-t330).*(-3.0./2.0);
t326 = t155.*t324;
t335 = t304+t332;
t338 = t318+t337;
t333 = t267+t289+t326;
et1 = t136.*((t23.*t63.*t94.*t167)./2.0+(t23.*t63.*t102.*t167)./2.0+t33.*t65.*t97.*t160.*(3.0./2.0)+(t35.*t66.*t97.*t161)./2.0+(t33.*t65.*t104.*t160)./2.0+(t35.*t66.*t104.*t161)./6.0-(t23.*t62.*t94.*(t229-t234))./2.0-(t23.*t62.*t102.*(t229-t234))./2.0+t23.*t63.*t101.*t160.*t167-t63.*t87.*t103.*t144.*t170.*2.0+(t63.*t87.*t97.*t160.*t170)./2.0+t63.*t87.*t104.*t160.*t170.*(3.0./2.0)-t23.*t62.*t101.*t160.*(t229-t234)+t21.*t23.*t63.*t101.*t144.*t170+(t23.*t33.*t65.*t94.*t144.*t167)./2.0-(t21.*t23.*t63.*t94.*t160.*t170)./2.0+(t23.*t33.*t65.*t102.*t144.*t167)./2.0-(t21.*t23.*t63.*t102.*t160.*t170)./2.0+t23.*t33.*t65.*t101.*t161.*t167);
et2 = -t135.*(t21.*t23.*t63.*t167+t23.*t33.*t65.*t160.*3.0+t23.*t35.*t66.*t161+t23.*t63.*t160.*t167+t63.*t88.*t160.*t170-t21.*t23.*t62.*(t229-t234)-t23.*t62.*t160.*(t229-t234)+t21.*t23.*t63.*t144.*t170+t23.*t33.*t65.*t161.*t167-t23.*t63.*t83.*t160.*t170+t21.*t23.*t33.*t65.*t144.*t167);
et3 = t37.*t43.*t50.*t55.*t60.*t335.*(t136.*(-t274-t277+t280+t281+t293+t4.*t38.*t74.*t78.*t87.*t101+(t4.*t38.*t74.*t78.*t87.*t176)./2.0-(t4.*t21.*t38.*t74.*t78.*t94.*t167)./2.0-(t4.*t21.*t38.*t74.*t78.*t102.*t167)./2.0+(t4.*t38.*t74.*t78.*t87.*t102.*t155)./2.0-t4.*t21.*t38.*t74.*t78.*t101.*t160.*t167)+t135.*(t249-t270-t4.*t38.*t74.*t78.*t87+t4.*t38.*t74.*t78.*t83.*t167-t4.*t38.*t74.*t78.*t87.*t163+t4.*t21.*t38.*t74.*t78.*t160.*t167)).*(4.5e+1./2.0)-t37.*t43.*t50.*t55.*t60.*(et1+et2).*(t291-t303).*(4.5e+1./2.0);
et4 = t37.*t38.*t43.*t50.*t55.*t90.*t160.*(t155.*(t279.*t309-t93.*(t113.*t283+t33.*t64.*t232.*8.0-t33.*t64.*t133.*(t119-t153).*6.0)-t86.*t143.*t168.*2.0+t82.*t143.*t166.*t279.*4.0)-t23.*t143.*t168.*8.0+t33.*t64.*t160.*t324+t143.*t160.*t166.*t238.*8.0-t33.*t64.*t89.*t93.*t132.*t142.*t156.*1.2e+1+(t35.*t89.*t93.*t132.*t142.*t155.*1.8e+1)./t40-t23.*t33.*t64.*t89.*t132.*t133.*t142.*t167.*6.0).*(4.5e+1./8.0)+t7.*t37.*t43.*t50.*t55.*t90.*t160.*t333.*(4.5e+1./4.0)-t37.*t43.*t50.*t55.*t61.*t62.*t161.*t333.*(4.5e+1./8.0);
et5 = t135.*(t23.*t63.*t160-t21.*t23.*t57.*t167+t23.*t33.*t65.*t161-t23.*t57.*t160.*t167+t23.*t63.*t161.*t167+t57.*t88.*t160.*t170-t21.*t23.*t56.*(t209-t233)-t23.*t56.*t160.*(t209-t233)+t21.*t23.*t57.*t144.*t170+t21.*t23.*t63.*t144.*t167-t23.*t57.*t83.*t160.*t170);
et6 = -t136.*((t63.*t97.*t160)./2.0+(t63.*t104.*t160)./6.0-(t23.*t57.*t94.*t167)./2.0-(t23.*t57.*t102.*t167)./2.0+(t33.*t65.*t97.*t161)./2.0+(t33.*t65.*t104.*t161)./6.0-(t23.*t56.*t94.*(t209-t233))./2.0-(t23.*t56.*t102.*(t209-t233))./2.0+(t23.*t63.*t94.*t144.*t167)./2.0+(t23.*t63.*t102.*t144.*t167)./2.0-t23.*t57.*t101.*t160.*t167+t23.*t63.*t101.*t161.*t167-t57.*t87.*t103.*t144.*t170.*2.0+(t57.*t87.*t97.*t160.*t170)./2.0+t57.*t87.*t104.*t160.*t170.*(3.0./2.0)-t23.*t56.*t101.*t160.*(t209-t233)+t21.*t23.*t57.*t101.*t144.*t170-(t21.*t23.*t57.*t94.*t160.*t170)./2.0-(t21.*t23.*t57.*t102.*t160.*t170)./2.0)-t36.*t59.*t197.*1.8e+1+t36.*t59.*t294.*1.8e+1;
et7 = t36.*t58.*t299.*-1.2e+1+t36.*t58.*t328.*1.2e+1;
et8 = (t6.*t42.*t49.*t53.*t59.*(t266-t297).*(9.0./2.0)+t6.*t42.*t49.*t53.*t58.*(t254-t296+t302-t330).*(3.0./2.0)).*(t298-t6.*t42.*t49.*t53.*t58.*(t136.*t327-t36.*t58.*t288.*6.0).*(3.0./2.0)).*1.0e+1-t37.*t43.*t50.*t55.*t59.^2.*t161.*t333.*(4.5e+1./8.0)+t37.*t38.*t43.*t50.*t55.*t57.^5.*t160.*t333.*(4.05e+2./8.0)-t6.*t42.*t49.*t53.*t58.*t136.*t288.*(t6.*t42.*t49.*t53.*t58.*(et5+et6+et7).*(-3.0./2.0)+t6.*t42.*t49.*t53.*t59.*(t254-t296+t302-t330).*9.0+t6.*t42.*t49.*t53.*t56.^5.*(t266-t297).*1.8e+1).*1.5e+1;
et9 = t37.*t38.*t43.*t50.*t55.*t90.*t160.*(t155.*(t278.*t309+t93.*(t171+t114.*(t36.*t58.*8.4e+1+t36.*t58.*t132.*8.4e+1)-t98.*t283+t133.*(t171-t107.*(t36.*t58.*1.6e+1-t36.*t58.*t132.*2.0e+1).*3.0+t4.*t63.*(t119-t153).*6.0)-t4.*t63.*t232.*8.0)-t86.*t143.*t165.*2.0+t82.*t157.*t166.*t278+t36.*t58.*t86.*t135.*t166.*2.4e+1)-t23.*t143.*t165.*8.0+t4.*t63.*t160.*t324+t143.*t160.*t166.*t237.*8.0+t23.*t36.*t58.*t135.*t166.*9.6e+1+t36.*t56.*t63.*t89.*t93.*t132.*t156.*1.8e+2+t36.*t56.*t63.*t89.*t93.*t142.*t156.*1.2e+1-t4.*t63.*t89.*t93.*t132.*t142.*t156.*1.2e+1+t34.*t65.*t89.*t93.*t132.*t142.*t155.*1.8e+1-t4.*t23.*t63.*t89.*t132.*t133.*t142.*t167.*6.0).*(4.5e+1./8.0);
et10 = ((t6.*t42.*t49.*t53.*t59.*(t266-t297).*(9.0./2.0)+t6.*t42.*t49.*t53.*t58.*(t254-t296+t302-t330).*(3.0./2.0)).*(-1.0./6.0))./t53-(t30.*(et8+et9))./7.2e+2;
et11 = (t10.*t40.*t75.*(-t155.*(-t138.*t315+t320.*(t4.*t22.*t73.*t77.*t160.*6.0+t4.*t22.*t73.*t77.*t81.*t160.*2.0)+t62.*t81.*t319+t4.*t62.*t81.*t316+t22.*t33.*t62.*t73.*t77.*t84.*t160.*t319.*2.0)-t22.*(t252.*4.0+t5.*t46.*t57.*t128-t5.*t46.*t57.*t148+t5.*t57.*t67.*t127-t5.*t57.*t67.*t149)+t20.*t62.*t319.*4.0+t22.*t137.*t315+t4.*t20.*t62.*t316.*4.0+t4.*t22.*t63.*t320.*8.0+t4.*t63.*t160.*(t321-t323)+t4.*t20.*t22.*t73.*t77.*t160.*t313+t33.*t62.*t73.*t77.*t85.*t160.*t319.*4.0-t4.*t20.*t22.*t73.*t77.*t137.*t160.*t320))./4.0;
et12 = t37.*t43.*t50.*t55.*t60.*(t135.*(t249-t270+t4.*t38.*t74.*t78.*3.0+t34.*t38.*t75.*t79.*t83.*2.0+t4.*t38.*t74.*t78.*t163.*3.0+t21.*t34.*t38.*t75.*t79.*t160.*2.0)-t136.*(t274+t277-t280-t281-t293+t4.*t38.*t74.*t78.*t101.*3.0+t4.*t38.*t74.*t78.*t176.*(3.0./2.0)+t21.*t34.*t38.*t75.*t79.*t94+t21.*t34.*t38.*t75.*t79.*t102+t4.*t38.*t74.*t78.*t102.*t155.*(3.0./2.0)+t21.*t34.*t38.*t75.*t79.*t101.*t160.*2.0)).*(t291-t303).*(4.5e+1./2.0);
et13 = t37.*t43.*t50.*t55.*t60.*t335.*(t135.*(t41.*t75.*t79.*t164.*2.0+t33.*t41.*t76.*t80.*t164-t21.*t23.*t41.*t75.*t79.*t133.*2.0)-t136.*(t33.*t41.*t76.*t80.*t103.*2.0+t33.*t41.*t76.*t80.*t198+t41.*t75.*t79.*t101.*t164.*2.0-t23.*t41.*t75.*t79.*t94.*t133-t23.*t41.*t75.*t79.*t102.*t133+t33.*t41.*t76.*t80.*t104.*t155.*(3.0./2.0))).*(4.5e+1./2.0)-t37.*t38.*t43.*t50.*t55.*t90.*t160.*(t155.*(t162.*t309+t82.*t152.*t157.*t166)+t21.*t143.*t152.*t166.*8.0+t89.*t107.*t132.*t142.*t156.*t162.*6.0).*(4.5e+1./8.0);
mt1 = [(t30.*(et3+et4))./7.2e+2-(t10.*t40.*t54.*t75.*(-t155.*(t320.*(t22.*t33.*t62.*t73.*t77.*t160.*6.0+t22.*t33.*t62.*t73.*t77.*t81.*t160.*2.0)+t4.*t63.*t81.*t319+t22.*t34.*t63.*t73.*t77.*t84.*t160.*t319.*2.0)+t4.*t20.*t63.*t319.*4.0+t22.*t33.*t64.*t320.*8.0+t33.*t64.*t160.*(t321-t323)+t34.*t63.*t73.*t77.*t85.*t160.*t319.*4.0+t20.*t22.*t33.*t62.*t73.*t77.*t160.*t313-t20.*t22.*t33.*t62.*t73.*t77.*t137.*t160.*t320).*3.0e+1+(t10.*t54.*t75.*(t22.*t313+t155.*(t321-t323)-t22.*t137.*t320-t4.*t20.*t62.*t319.*4.0).*1.5e+2)./t65)./(t54.*1.2e+2)+(t6.*t42.*t49.*t58.*t335)./4.0;et10+et11];
mt2 = [t10.*t40.*t75.*(t155.*(t138.*t311-t4.*t62.*t81.*t312)-t22.*(t243.*4.0+t46.*t56.*t128-t46.*t56.*t148+t56.*t67.*t127-t56.*t67.*t149)+t22.*t137.*t311+t4.*t20.*t62.*t312.*4.0).*(-1.0./4.0)+(t6.*t42.*t49.*t58.*(t5.*t57.*t197.*6.0-t5.*t57.*t294.*6.0))./4.0+(t37.*t38.*t43.*t50.*t90.*t160.*(t155.*(t93.*(t169+t114.*(t5.*t57.*8.4e+1+t5.*t57.*t132.*8.4e+1)+t133.*(t169-t107.*(t5.*t57.*1.6e+1-t5.*t57.*t132.*2.0e+1).*3.0))+t5.*t57.*t86.*t135.*t166.*2.4e+1)+t5.*t23.*t57.*t135.*t166.*9.6e+1+t5.*t63.*t89.*t93.*t132.*t156.*1.8e+2+t5.*t63.*t89.*t93.*t142.*t156.*1.2e+1))./1.28e+2];
mt3 = [(t30.*(et12+et13))./7.2e+2-(t10.*t40.*t75.*(t155.*(t320.*(t38.*t73.*t77.*6.0+t38.*t73.*t77.*t81.*2.0-6.0)+t4.*t7.*t73.*t77.*t84.*t319.*2.0)-t20.*t38.*t73.*t77.*t313-t4.*t7.*t22.*t73.*t77.*t319.*4.0+t20.*t38.*t73.*t77.*t137.*t320))./4.0-(t6.*t42.*t49.*t58.*(t291-t303))./4.0];
mt4 = [t30.*(t6.*t42.*t49.*t53.*t58.*t136.*(t6.*t42.*t49.*t53.*t59.*(t266-t297).*(9.0./2.0)+t6.*t42.*t49.*t53.*t58.*(t254-t296+t302-t330).*(3.0./2.0)).*(t103.*2.0+t97.*t155.*2.0+t104.*t155.*(2.0./3.0)).*-1.5e+1+t6.*t42.*t49.*t53.*t58.*t136.*t288.*(t298+t6.*t36.*t42.*t49.*t53.*t60.*t288.*9.0-t6.*t42.*t49.*t53.*t58.*t136.*t327.*(3.0./2.0)).*1.5e+1+t37.*t43.*t50.*t55.*t61.*t93.*t132.*t133.*t142.*sin(t47).*(1.35e+2./2.0)).*(-1.0./7.2e+2)+(t10.*t40.*t75.*(t22.*(t173.*4.0+t192+t193-t227+t228)-t155.*(t138.*t319+t4.*t62.*t81.*t320)-t22.*t137.*t319+t4.*t20.*t62.*t320.*4.0))./4.0+(t6.*t42.*t49.*t58.*t136.*t288)./4.0];
mt5 = [(t10.*t40.*t75.*(t22.*(t192-t193+t227+t228)-t155.*(t138.*(t217-t218+t240+t245)+t4.*t62.*t81.*(t216+t222+t235-t244))-t22.*t137.*(t217-t218+t240+t245)+t4.*t20.*t62.*(t216+t222+t235-t244).*4.0))./4.0;0.0;t10.*t40.*t75.*(t155.*(t138.*t314-t4.*t62.*t81.*t317)-t22.*(t174.*4.0+t192+t193+t212+t227)+t22.*t137.*t314+t4.*t20.*t62.*t317.*4.0).*(-1.0./4.0)];
getFinalDelOffset = [mt1;mt2;mt3;mt4;mt5];
