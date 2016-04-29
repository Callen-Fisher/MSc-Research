function F = FfunctionEKF(dth1,dth2,dth3,dth4,th1,th2,th3,th4)
%FFUNCTIONEKF
%    F = FFUNCTIONEKF(DTH1,DTH2,DTH3,DTH4,TH1,TH2,TH3,TH4)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    14-Jun-2015 17:58:57

t2 = th1.*2.0;
t3 = th4.*2.0;
t13 = th1-th2;
t4 = sin(t13);
t16 = th1-th3;
t5 = sin(t16);
t19 = th1-th4;
t6 = sin(t19);
t22 = th2-th3;
t7 = sin(t22);
t25 = th2-th4;
t8 = sin(t25);
t28 = th3-th4;
t9 = sin(t28);
t31 = th1-th2+th3-th4;
t10 = sin(t31);
t33 = th1-th2-th3+th4;
t11 = sin(t33);
t12 = th3.*2.0;
t14 = t4.^2;
t15 = t14.*1.413644886e9;
t17 = t5.^2;
t18 = t17.*2.911970817e10;
t20 = t6.^2;
t21 = t20.*1.17322695e9;
t23 = t7.^2;
t24 = t23.*1.081443069e10;
t26 = t8.^2;
t27 = t26.*4.3571115e8;
t29 = t9.^2;
t30 = t29.*9.8250151959e11;
t32 = t10.^2;
t34 = t11.^2;
t37 = t32.*3.5724645e7;
t38 = t34.*3.5724645e7;
t35 = t15+t18+t21+t24+t27+t30-t37-t38+6.924711584315e12;
t36 = 1.0./t35;
t39 = dth1.^2;
t41 = th2.*2.0;
t40 = t2-t3+t12-t41;
t42 = t2-t12;
t43 = t2-t3;
t44 = dth2.^2;
t45 = -t12+th1+th2;
t46 = -t3+th1+th2;
t47 = dth3.^2;
t48 = -t3+th1+th3;
t49 = dth4.^2;
t50 = -t12+th1+th4;
t51 = t3-t12+th1-th2;
t52 = -t3+t12+th1-th2;
t53 = cos(t13);
t54 = cos(t16);
t55 = cos(t19);
t56 = t2+t3-t12-t41;
t57 = sin(t40);
t58 = t2-t41;
t59 = sin(t42);
t60 = sin(t43);
t61 = sin(t45);
t62 = sin(t46);
t63 = -t41+th1+th3;
t64 = sin(t48);
t65 = -t41+th1+th4;
t66 = sin(t50);
t67 = sin(t51);
t68 = sin(t52);
t69 = t3-t41+th1-th3;
t70 = t12-t41+th1-th4;
t71 = t44.*t53.*2.7950393321e11;
t72 = cos(t56);
t73 = cos(t40);
t74 = cos(t58);
t75 = t39.*t74.*3.92679135e9;
t76 = cos(t45);
t77 = t44.*t76.*2.4266423475e10;
t78 = cos(t46);
t79 = t44.*t78.*9.77689125e8;
t80 = cos(t63);
t81 = cos(t65);
t82 = cos(t51);
t83 = cos(t52);
t84 = cos(t69);
t85 = cos(t70);
t86 = cos(t31);
t87 = cos(t33);
t88 = t4.*t53.*2.827289772e9;
t89 = 1.0./t35.^2;
t90 = t4.*t44.*2.7950393321e11;
t91 = t5.*t47.*7.2757812987e10;
t92 = t6.*t49.*1.24186747525e10;
t93 = sin(t56);
t94 = sin(t58);
t95 = t39.*t94.*1.963395675e9;
t96 = t39.*t59.*4.0444039125e10;
t97 = t39.*t60.*1.629481875e9;
t98 = t44.*t61.*2.4266423475e10;
t99 = t44.*t62.*9.77689125e8;
t100 = sin(t63);
t101 = t47.*t100.*2.83242987e8;
t102 = t47.*t64.*1.177597305e9;
t103 = sin(t65);
t104 = t49.*t103.*4.97612025e7;
t105 = t49.*t66.*8.1189392175e9;
t106 = sin(t69);
t107 = t47.*t106.*4.370355e6;
t108 = sin(t70);
t109 = t49.*t108.*3.01313925e7;
t125 = t39.*t93.*4.96175625e7;
t126 = t39.*t57.*4.96175625e7;
t127 = t44.*t67.*8.05147695e9;
t128 = t44.*t68.*8.05147695e9;
t110 = t90+t91+t92+t95+t96+t97+t98+t99+t101+t102+t104+t105+t107+t109-t125-t126-t127-t128;
t111 = t47.*t54.*7.2757812987e10;
t112 = cos(t42);
t113 = t39.*t112.*8.088807825e10;
t114 = t47.*t80.*2.83242987e8;
t115 = cos(t48);
t116 = t47.*t115.*1.177597305e9;
t117 = cos(t50);
t118 = t47.*t84.*4.370355e6;
t119 = t49.*t85.*6.0262785e7;
t120 = t10.*t86.*7.144929e7;
t121 = t11.*t87.*7.144929e7;
t122 = t5.*t54.*5.823941634e10;
t123 = cos(t22);
t124 = t7.*t123.*2.162886138e10;
t129 = t49.*t55.*1.24186747525e10;
t130 = t39.*t73.*9.9235125e7;
t131 = cos(t43);
t132 = t39.*t131.*3.25896375e9;
t133 = t49.*t81.*4.97612025e7;
t134 = t49.*t117.*8.1189392175e9;
t135 = t44.*t83.*1.61029539e10;
t136 = t47.*t84.*8.74071e6;
t137 = t49.*t85.*3.01313925e7;
t138 = t6.*t55.*2.3464539e9;
t139 = cos(t25);
t140 = t8.*t139.*8.714223e8;
t141 = cos(t28);
t142 = t2-th2-th3;
t143 = t2-th2-th4;
t144 = t2-t3-th2+th3;
t145 = t2-t12-th2+th4;
t146 = sin(t142);
t147 = sin(t143);
t148 = -t12+t41;
t149 = sin(t148);
t150 = -t3+t41;
t151 = sin(t150);
t152 = -t3+th2+th3;
t153 = sin(t152);
t154 = -t12+th2+th4;
t155 = sin(t154);
t156 = sin(t144);
t157 = sin(t145);
t158 = t39.*t53.*5.849044645e10;
t159 = cos(t142);
t160 = cos(t143);
t161 = t44.*t74.*7.8535827e8;
t162 = cos(t144);
t163 = cos(t145);
t164 = -t88+t120+t121+t124+t140;
t165 = t7.*t47.*9.066244329e9;
t166 = t8.*t49.*1.5477695175e9;
t167 = t44.*t93.*9.9235125e6;
t168 = t44.*t57.*9.9235125e6;
t169 = t44.*t149.*3.004008525e9;
t170 = t44.*t151.*1.21030875e8;
t171 = t39.*t61.*5.006680875e9;
t172 = t39.*t62.*2.01718125e8;
t173 = t47.*t153.*1.46693835e8;
t174 = t49.*t155.*1.0113799725e9;
t175 = t39.*t67.*1.68198525e9;
t176 = t39.*t68.*1.68198525e9;
t188 = t4.*t39.*5.849044645e10;
t189 = t47.*t146.*9.4414329e7;
t190 = t49.*t147.*1.65870675e7;
t191 = t44.*t94.*3.92679135e8;
t192 = t47.*t156.*1.456785e6;
t193 = t49.*t157.*1.00437975e7;
t177 = t165+t166+t167+t168+t169+t170+t171+t172+t173+t174+t175+t176-t188-t189-t190-t191-t192-t193;
t178 = t47.*t123.*9.066244329e9;
t179 = t47.*t159.*9.4414329e7;
t180 = cos(t148);
t181 = t44.*t180.*6.00801705e9;
t182 = cos(t152);
t183 = t47.*t182.*1.46693835e8;
t184 = cos(t154);
t185 = t47.*t162.*1.456785e6;
t186 = t49.*t163.*2.0087595e7;
t187 = t9.*t141.*1.96500303918e12;
t194 = t49.*t139.*1.5477695175e9;
t195 = t44.*t72.*1.9847025e7;
t196 = t49.*t160.*1.65870675e7;
t197 = cos(t150);
t198 = t44.*t197.*2.4206175e8;
t199 = t49.*t184.*1.0113799725e9;
t200 = t39.*t82.*3.3639705e9;
t201 = t47.*t162.*2.91357e6;
t202 = t49.*t163.*1.00437975e7;
t203 = -t120+t121+t138+t140+t187;
t204 = t14.*2.827289772e9;
t205 = t17.*5.823941634e10;
t206 = t20.*2.3464539e9;
t207 = t23.*2.162886138e10;
t208 = t26.*8.714223e8;
t209 = t29.*1.96500303918e12;
t212 = t32.*7.144929e7;
t213 = t34.*7.144929e7;
t210 = t204+t205+t206+t207+t208+t209-t212-t213+1.384942316863e13;
t211 = 1.0./t210;
t214 = t2-th3-th4;
t215 = t2-t41+th3-th4;
t216 = t2-t41-th3+th4;
t217 = sin(t214);
t218 = t41-th3-th4;
t219 = sin(t218);
t220 = -t3+t12;
t221 = sin(t220);
t222 = sin(t215);
t223 = sin(t216);
t224 = t47.*t73.*2.64627e6;
t225 = cos(t215);
t226 = t49.*t225.*1.188594e7;
t227 = cos(t216);
t228 = t4.*t53.*5.654579544e9;
t229 = 1.0./t210.^2;
t230 = t5.*t39.*6.14505700375e11;
t231 = t7.*t44.*3.65880787125e11;
t232 = t47.*t57.*1.323135e6;
t233 = t44.*t146.*4.490194275e9;
t234 = t49.*t217.*2.99585475e8;
t235 = t49.*t219.*1.11259575e8;
t236 = t47.*t59.*1.07850771e9;
t237 = t47.*t149.*4.0053447e8;
t238 = t39.*t106.*1.83211875e8;
t239 = t49.*t222.*5.94297e6;
t254 = t9.*t49.*1.3316533994e11;
t255 = t47.*t93.*1.323135e6;
t256 = t47.*t221.*3.638894517e10;
t257 = t39.*t100.*2.779268625e9;
t258 = t39.*t64.*5.0353314375e10;
t259 = t44.*t153.*3.0025918125e10;
t260 = t44.*t156.*2.95997625e8;
t261 = t49.*t223.*5.94297e6;
t240 = t230+t231+t232+t233+t234+t235+t236+t237+t238+t239-t254-t255-t256-t257-t258-t259-t260-t261;
t241 = t39.*t54.*6.14505700375e11;
t242 = t44.*t123.*3.65880787125e11;
t243 = t47.*t72.*2.64627e6;
t244 = cos(t214);
t245 = cos(t218);
t246 = t47.*t112.*2.15701542e9;
t247 = t47.*t180.*8.0106894e8;
t248 = t39.*t84.*1.83211875e8;
t249 = t44.*t162.*2.95997625e8;
t250 = t10.*t86.*1.4289858e8;
t251 = t11.*t87.*1.4289858e8;
t252 = t5.*t54.*1.1647883268e11;
t253 = t7.*t123.*4.325772276e10;
t262 = t49.*t141.*1.3316533994e11;
t263 = t49.*t244.*2.99585475e8;
t264 = t49.*t245.*1.11259575e8;
t265 = cos(t220);
t266 = t47.*t265.*7.277789034e10;
t267 = t6.*t55.*4.6929078e9;
t268 = t8.*t139.*1.7428446e9;
t269 = t88-t120-t121+t122+t138;
t270 = t49.*t72.*5.9541075e6;
t271 = t47.*t227.*5.5106568e7;
t272 = t6.*t39.*5.919714574375e11;
t273 = t8.*t44.*3.528152683125e11;
t274 = t9.*t47.*7.22199382848e11;
t275 = t49.*t93.*2.97705375e6;
t276 = t44.*t147.*3.7663455375e9;
t277 = t47.*t217.*3.51968085e8;
t278 = t47.*t219.*1.30713345e8;
t279 = t49.*t60.*9.77689125e7;
t280 = t49.*t151.*3.63092625e7;
t281 = t49.*t221.*8.18751266325e10;
t282 = t39.*t108.*1.4840161875e9;
t283 = t47.*t223.*2.7553284e7;
t287 = t49.*t57.*2.97705375e6;
t288 = t39.*t103.*2.3312323125e9;
t289 = t39.*t66.*4.078618464375e11;
t290 = t44.*t155.*2.432099368125e11;
t291 = t44.*t157.*2.3975807625e9;
t292 = t47.*t222.*2.7553284e7;
t284 = t272+t273+t274+t275+t276+t277+t278+t279+t280+t281+t282+t283-t287-t288-t289-t290-t291-t292;
t285 = t49.*t73.*5.9541075e6;
t286 = t120-t121+t122+t124-t187;
t293 = t39.*t55.*5.919714574375e11;
t294 = t44.*t139.*3.528152683125e11;
t295 = t47.*t244.*3.51968085e8;
t296 = t47.*t245.*1.30713345e8;
t297 = t49.*t131.*1.95537825e8;
t298 = t49.*t197.*7.2618525e7;
t299 = t39.*t85.*1.4840161875e9;
t300 = t44.*t163.*2.3975807625e9;
t301 = t47.*t225.*2.7553284e7;
t302 = t47.*t227.*2.7553284e7;
F = reshape([t36.*(dth1.*sin(t2+t3-th2.*2.0-th3.*2.0).*-9.9235125e7-dth1.*t57.*9.9235125e7+dth1.*t59.*8.088807825e10+dth1.*t60.*3.25896375e9+dth1.*sin(t2-th2.*2.0).*3.92679135e9).*(-9.0./2.5e1),t36.*(dth1.*t4.*-1.169808929e11+dth1.*t61.*1.001336175e10+dth1.*t62.*4.0343625e8+dth1.*t67.*3.3639705e9+dth1.*t68.*3.3639705e9).*(-9.0./5.0),t211.*(dth1.*t5.*1.22901140075e12-dth1.*t64.*1.0070662875e11-dth1.*t100.*5.55853725e9+dth1.*t106.*3.6642375e8).*2.7e1,t36.*(dth1.*t6.*1.183942914875e12-dth1.*t66.*8.15723692875e11-dth1.*t103.*4.662464625e9+dth1.*t108.*2.968032375e9).*6.0,1.0,0.0,0.0,0.0,t36.*(dth2.*t4.*5.5900786642e11+dth2.*t61.*4.853284695e10+dth2.*t62.*1.95537825e9-dth2.*t67.*1.61029539e10-dth2.*t68.*1.61029539e10).*(-9.0./2.5e1),t36.*(dth2.*t57.*1.9847025e7+dth2.*t93.*1.9847025e7-dth2.*t94.*7.8535827e8+dth2.*t149.*6.00801705e9+dth2.*t151.*2.4206175e8).*(-9.0./5.0),t211.*(dth2.*t7.*7.3176157425e11+dth2.*t146.*8.98038855e9-dth2.*t153.*6.005183625e10-dth2.*t156.*5.9199525e8).*2.7e1,t36.*(dth2.*t8.*7.05630536625e11+dth2.*t147.*7.532691075e9-dth2.*t155.*4.86419873625e11-dth2.*t157.*4.795161525e9).*6.0,0.0,1.0,0.0,0.0,t36.*(dth3.*sin(t3+th1-th2.*2.0-th3).*8.74071e6+dth3.*t5.*1.45515625974e11+dth3.*t64.*2.35519461e9+dth3.*sin(th1-th2.*2.0+th3).*5.66485974e8).*(-9.0./2.5e1),t36.*(dth3.*t7.*1.8132488658e10-dth3.*t146.*1.88828658e8+dth3.*t153.*2.9338767e8-dth3.*t156.*2.91357e6).*(-9.0./5.0),t211.*(dth3.*t57.*2.64627e6+dth3.*t59.*2.15701542e9-dth3.*t93.*2.64627e6+dth3.*t149.*8.0106894e8-dth3.*t221.*7.277789034e10).*2.7e1,t36.*(dth3.*t9.*1.444398765696e12+dth3.*t217.*7.0393617e8+dth3.*t219.*2.6142669e8-dth3.*t222.*5.5106568e7+dth3.*t223.*5.5106568e7).*6.0,0.0,0.0,1.0,0.0,t36.*(dth4.*sin(t12+th1-th2.*2.0-th4).*6.0262785e7+dth4.*t6.*2.4837349505e10+dth4.*t66.*1.6237878435e10+dth4.*sin(th1-th2.*2.0+th4).*9.9522405e7).*(-9.0./2.5e1),t36.*(dth4.*t8.*3.095539035e9-dth4.*t147.*3.3174135e7+dth4.*t155.*2.022759945e9-dth4.*t157.*2.0087595e7).*(-9.0./5.0),t211.*(dth4.*t9.*-2.6633067988e11+dth4.*t217.*5.9917095e8+dth4.*t219.*2.2251915e8+dth4.*t222.*1.188594e7-dth4.*t223.*1.188594e7).*2.7e1,t36.*(dth4.*t57.*(-5.9541075e6)+dth4.*t60.*1.95537825e8+dth4.*t93.*5.9541075e6+dth4.*t151.*7.2618525e7+dth4.*t221.*1.63750253265e11).*6.0,0.0,0.0,0.0,1.0,t36.*(t71+t75+t77+t79+t111+t113+t114+t116+t118+t129+t132+t133+t134+t137-t39.*t72.*9.9235125e7-t39.*t73.*9.9235125e7-t44.*t82.*8.05147695e9-t44.*t83.*8.05147695e9).*(-9.0./2.5e1)+t89.*t110.*(t88+t122+t138-t10.*t86.*7.144929e7-t11.*t87.*7.144929e7).*(9.0./2.5e1),t36.*(t158+t161+t186+t201-t39.*t76.*5.006680875e9-t44.*t72.*1.9847025e7-t39.*t78.*2.01718125e8-t44.*t73.*1.9847025e7-t39.*t82.*1.68198525e9-t39.*t83.*1.68198525e9+t47.*t159.*1.88828658e8+t49.*t160.*3.3174135e7).*(9.0./5.0)+t89.*t177.*t269.*(9.0./5.0),t211.*(t224+t226+t241+t246+t248-t39.*t80.*2.779268625e9-t47.*t72.*2.64627e6-t39.*t115.*5.0353314375e10+t44.*t159.*8.98038855e9-t44.*t162.*5.9199525e8-t49.*t227.*1.188594e7+t49.*t244.*5.9917095e8).*2.7e1-t229.*t240.*(t228+t252+t267-t10.*t86.*1.4289858e8-t11.*t87.*1.4289858e8).*2.7e1,t36.*(t270+t271+t293+t297+t299-t39.*t81.*2.3312323125e9-t49.*t73.*5.9541075e6-t39.*t117.*4.078618464375e11+t44.*t160.*7.532691075e9-t44.*t163.*4.795161525e9-t47.*t225.*5.5106568e7+t47.*t244.*7.0393617e8).*6.0-t89.*t269.*t284.*6.0,0.0,0.0,0.0,0.0,t36.*(t71+t75-t77-t79+t119+t136-t39.*t72.*9.9235125e7-t39.*t73.*9.9235125e7-t44.*t82.*8.05147695e9-t44.*t83.*8.05147695e9+t47.*t80.*5.66485974e8+t49.*t81.*9.9522405e7).*(9.0./2.5e1)+t89.*t110.*t164.*(9.0./2.5e1),t36.*(t158+t161+t178+t179+t181+t183+t185+t194+t196+t198+t199+t202+t39.*t76.*5.006680875e9-t44.*t72.*1.9847025e7+t39.*t78.*2.01718125e8-t44.*t73.*1.9847025e7-t39.*t82.*1.68198525e9-t39.*t83.*1.68198525e9).*(-9.0./5.0)+t89.*t164.*t177.*(9.0./5.0),t211.*(-t224-t226+t242+t243+t247+t249+t39.*t80.*5.55853725e9-t39.*t84.*3.6642375e8-t44.*t159.*4.490194275e9-t44.*t182.*3.0025918125e10+t49.*t227.*1.188594e7+t49.*t245.*2.2251915e8).*2.7e1-t229.*t240.*(-t228+t250+t251+t253+t268).*2.7e1,t36.*(-t270-t271+t285+t294+t298+t300+t39.*t81.*4.662464625e9-t39.*t85.*2.968032375e9-t44.*t160.*3.7663455375e9-t44.*t184.*2.432099368125e11+t47.*t225.*5.5106568e7+t47.*t245.*2.6142669e8).*6.0-t89.*t164.*t284.*6.0,0.0,0.0,0.0,0.0,t36.*(t111+t113-t114-t116+t118-t119+t130+t135-t39.*t72.*9.9235125e7+t44.*t76.*4.853284695e10-t44.*t82.*1.61029539e10+t49.*t117.*1.6237878435e10).*(9.0./2.5e1)-t89.*t110.*(t120-t121+t122+t124-t9.*t141.*1.96500303918e12).*(9.0./2.5e1),t36.*(t178-t179+t181-t183+t185-t186+t195+t200+t39.*t76.*1.001336175e10-t44.*t73.*1.9847025e7-t39.*t83.*3.3639705e9+t49.*t184.*2.022759945e9).*(9.0./5.0)-t89.*t177.*t286.*(9.0./5.0),t211.*(-t224+t241+t242-t243+t246+t247+t248+t249+t262+t263+t264+t266+t39.*t80.*2.779268625e9+t39.*t115.*5.0353314375e10+t44.*t159.*4.490194275e9+t44.*t182.*3.0025918125e10-t49.*t225.*5.94297e6-t49.*t227.*5.94297e6).*-2.7e1+t229.*t240.*(t250-t251+t252+t253-t9.*t141.*3.93000607836e12).*2.7e1,t36.*(t270+t285+t295+t296+t301+t302-t39.*t85.*2.968032375e9-t39.*t117.*8.15723692875e11-t47.*t141.*7.22199382848e11-t44.*t163.*4.795161525e9-t44.*t184.*4.86419873625e11-t49.*t265.*1.63750253265e11).*-6.0+t89.*t284.*t286.*6.0,0.0,0.0,0.0,0.0,t36.*(t129-t130+t132-t133-t134-t135-t136+t137+t39.*t72.*9.9235125e7+t44.*t78.*1.95537825e9+t44.*t82.*1.61029539e10+t47.*t115.*2.35519461e9).*(9.0./2.5e1)-t89.*t110.*t203.*(9.0./2.5e1),t36.*(t194-t195-t196+t198-t199-t200-t201+t202+t39.*t78.*4.0343625e8+t44.*t73.*1.9847025e7+t39.*t83.*3.3639705e9+t47.*t182.*2.9338767e8).*(9.0./5.0)-t89.*t177.*t203.*(9.0./5.0),t211.*(t224+t243-t262+t263+t264-t266-t39.*t84.*3.6642375e8-t39.*t115.*1.0070662875e11-t44.*t162.*5.9199525e8-t44.*t182.*6.005183625e10+t49.*t225.*5.94297e6+t49.*t227.*5.94297e6).*-2.7e1+t229.*t240.*(-t250+t251+t267+t268+t9.*t141.*3.93000607836e12).*2.7e1,t36.*(-t270-t285+t293+t294+t295+t296+t297+t298+t299+t300-t301-t302+t39.*t81.*2.3312323125e9+t39.*t117.*4.078618464375e11+t47.*t141.*7.22199382848e11+t44.*t160.*3.7663455375e9+t44.*t184.*2.432099368125e11+t49.*t265.*1.63750253265e11).*-6.0+t89.*t203.*t284.*6.0,0.0,0.0,0.0,0.0],[8, 8]);
