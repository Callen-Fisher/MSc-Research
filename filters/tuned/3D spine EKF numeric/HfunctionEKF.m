function H = HfunctionEKF(ddph1,ddph2,ddph3,ddph4,ddth1,ddth2,ddth3,ddth4,dph1,dph2,dph3,dph4,dsi1,dsi2,dth1,dth2,dth3,dth4,ph1,ph2,ph3,ph4,si1,si2,th1,th2,th3,th4)
%HFUNCTIONEKF
%    H = HFUNCTIONEKF(DDPH1,DDPH2,DDPH3,DDPH4,DDTH1,DDTH2,DDTH3,DDTH4,DPH1,DPH2,DPH3,DPH4,DSI1,DSI2,DTH1,DTH2,DTH3,DTH4,PH1,PH2,PH3,PH4,SI1,SI2,TH1,TH2,TH3,TH4)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    13-Jun-2015 13:01:31

t2 = cos(th1);
t3 = cos(ph1);
t4 = sin(th1);
t5 = sin(ph1);
t6 = sqrt(2.0);
t7 = sin(si1);
t8 = cos(si1);
t9 = t3.*t8;
t10 = t4.*t5.*t7;
t11 = t9+t10;
t12 = t5.*t8;
t13 = t12-t3.*t4.*t7;
t14 = t5.*t7;
t15 = t3.*t4.*t8;
t16 = t14+t15;
t17 = t6.*t16.*(1.0./2.0);
t18 = t3.*t7;
t19 = t18-t4.*t5.*t8;
t20 = cos(th2);
t21 = cos(ph2);
t22 = sin(ph2);
t23 = sin(th2);
t24 = dph1.*t2.*t5.*(1.0./2.0);
t25 = dth1.*t3.*t4.*(1.0./2.0);
t26 = t24+t25;
t27 = dph1.*t3.*t4.*(1.0./2.0);
t28 = dth1.*t2.*t5.*(1.0./2.0);
t29 = t27+t28;
t30 = dph1.*t2.*t3.*(1.0./2.0);
t36 = dth1.*t4.*t5.*(1.0./2.0);
t31 = t30-t36;
t32 = dph1.*t4.*t5.*(1.0./2.0);
t38 = dth1.*t2.*t3.*(1.0./2.0);
t33 = t32-t38;
t34 = dth1.^2;
t35 = dth2.^2;
t37 = dph1.*t31;
t39 = ddph1.*t2.*t5.*(1.0./2.0);
t40 = ddth1.*t3.*t4.*(1.0./2.0);
t41 = dph1.*t26;
t42 = dth1.*t29;
t43 = ddth1.*t4.*t5.*(1.0./2.0);
t44 = dph2.*t21.*t23.*(3.0./1.0e1);
t45 = dth2.*t20.*t22.*(3.0./1.0e1);
t46 = t44+t45;
t47 = dph2.*t20.*t22.*(3.0./1.0e1);
t48 = dth2.*t21.*t23.*(3.0./1.0e1);
t49 = t47+t48;
t50 = dph2.*t22.*t23.*(3.0./1.0e1);
t56 = dth2.*t20.*t21.*(3.0./1.0e1);
t51 = t50-t56;
t52 = dph2.*t20.*t21.*(3.0./1.0e1);
t54 = dth2.*t22.*t23.*(3.0./1.0e1);
t53 = t52-t54;
t55 = dph2.*t53;
t57 = ddph2.*t20.*t22.*(3.0./1.0e1);
t58 = ddth2.*t21.*t23.*(3.0./1.0e1);
t63 = dth2.*t51;
t103 = dth1.*t33;
t59 = t37+t39+t40+t55+t57+t58-t63-t103;
t60 = dph2.*t49;
t61 = dth2.*t46;
t62 = ddth2.*t22.*t23.*(3.0./1.0e1);
t102 = ddph1.*t2.*t3.*(1.0./2.0);
t104 = ddph2.*t20.*t21.*(3.0./1.0e1);
t64 = t41+t42+t43+t60+t61+t62-t102-t104;
t65 = dph1.*t3.*t4;
t66 = dth1.*t2.*t5;
t67 = t65+t66;
t68 = dph1.*t4.*t5;
t69 = sin(si2);
t70 = cos(si2);
t71 = dph1.*t2.*t5;
t72 = dth1.*t3.*t4;
t73 = t71+t72;
t74 = t21.*t69;
t85 = t22.*t23.*t70;
t75 = t74-t85;
t76 = dph1.*t2.*t3;
t128 = dth1.*t4.*t5;
t77 = t76-t128;
t78 = t22.*t69;
t79 = t21.*t23.*t70;
t80 = t78+t79;
t81 = dph2.*t22.*t23.*(3.0./5.0);
t82 = dph2.*t21.*t23.*(3.0./5.0);
t83 = dth2.*t20.*t22.*(3.0./5.0);
t84 = t82+t83;
t86 = dph2.*t20.*t21.*(3.0./5.0);
t131 = dth2.*t22.*t23.*(3.0./5.0);
t87 = t86-t131;
t88 = dph2.*t20.*t22.*(3.0./5.0);
t89 = dth2.*t21.*t23.*(3.0./5.0);
t90 = t88+t89;
t91 = dph1.*t29;
t92 = dth1.*t26;
t93 = ddph1.*t4.*t5.*(1.0./2.0);
t132 = ddth1.*t2.*t3.*(1.0./2.0);
t94 = t91+t92+t93-t132;
t95 = dth1.*t31;
t96 = ddph1.*t3.*t4.*(1.0./2.0);
t97 = ddth1.*t2.*t5.*(1.0./2.0);
t133 = dph1.*t33;
t98 = t95+t96+t97-t133;
t99 = ddth1.*t4.*(1.0./2.0);
t100 = t2.*t34.*(1.0./2.0);
t101 = t99+t100;
t105 = ddth1.*t2.*(1.0./2.0);
t106 = ddth2.*t20.*(3.0./1.0e1);
t107 = dph2.*t46;
t108 = dth2.*t49;
t109 = ddph2.*t22.*t23.*(3.0./1.0e1);
t136 = ddth2.*t20.*t21.*(3.0./1.0e1);
t110 = t107+t108+t109-t136;
t111 = dth2.*t53;
t112 = ddph2.*t21.*t23.*(3.0./1.0e1);
t113 = ddth2.*t20.*t22.*(3.0./1.0e1);
t137 = dph2.*t51;
t114 = t111+t112+t113-t137;
t115 = ddth2.*t23.*(3.0./1.0e1);
t116 = t20.*t35.*(3.0./1.0e1);
t117 = t115+t116;
t118 = t55+t57+t58-t63;
t119 = t21.*t70;
t120 = t22.*t23.*t69;
t121 = t119+t120;
t122 = t22.*t70;
t124 = t21.*t23.*t69;
t123 = t122-t124;
t125 = t6.*t80.*(1.0./2.0);
t126 = cos(th3);
t127 = cos(ph3);
t129 = sin(ph3);
t130 = sin(th3);
t134 = t41+t42+t43-t102;
t135 = t37+t39+t40-t103;
t138 = t60+t61+t62-t104;
t139 = dth3.^2;
t140 = dph3.*t127.*t130.*(9.0./2.5e1);
t141 = dth3.*t126.*t129.*(9.0./2.5e1);
t142 = t140+t141;
t143 = dph3.*t126.*t129.*(9.0./2.5e1);
t144 = dth3.*t127.*t130.*(9.0./2.5e1);
t145 = t143+t144;
t146 = dph3.*t129.*t130.*(9.0./2.5e1);
t156 = dth3.*t126.*t127.*(9.0./2.5e1);
t147 = t146-t156;
t148 = dph3.*t126.*t127.*(9.0./2.5e1);
t154 = dth3.*t129.*t130.*(9.0./2.5e1);
t149 = t148-t154;
t150 = dph3.*t145;
t151 = dth3.*t142;
t152 = ddth3.*t129.*t130.*(9.0./2.5e1);
t160 = ddph3.*t126.*t127.*(9.0./2.5e1);
t153 = t41+t42+t43+t60+t61+t62-t102-t104+t150+t151+t152-t160;
t155 = dph3.*t149;
t157 = ddph3.*t126.*t129.*(9.0./2.5e1);
t158 = ddth3.*t127.*t130.*(9.0./2.5e1);
t161 = dth3.*t147;
t159 = t37+t39+t40+t55+t57+t58-t63-t103+t155+t157+t158-t161;
t162 = dph3.*t129.*t130.*(1.8e1./2.5e1);
t163 = dph3.*t127.*t130.*(1.8e1./2.5e1);
t164 = dth3.*t126.*t129.*(1.8e1./2.5e1);
t165 = t163+t164;
t166 = dph3.*t126.*t127.*(1.8e1./2.5e1);
t190 = dth3.*t129.*t130.*(1.8e1./2.5e1);
t167 = t166-t190;
t168 = dph3.*t126.*t129.*(1.8e1./2.5e1);
t169 = dth3.*t127.*t130.*(1.8e1./2.5e1);
t170 = t168+t169;
t171 = ddth3.*t130.*(9.0./2.5e1);
t172 = t126.*t139.*(9.0./2.5e1);
t173 = t171+t172;
t174 = ddth3.*t126.*(9.0./2.5e1);
t175 = dph3.*t142;
t176 = dth3.*t145;
t177 = ddph3.*t129.*t130.*(9.0./2.5e1);
t191 = ddth3.*t126.*t127.*(9.0./2.5e1);
t178 = t175+t176+t177-t191;
t179 = dth3.*t149;
t180 = ddph3.*t127.*t130.*(9.0./2.5e1);
t181 = ddth3.*t126.*t129.*(9.0./2.5e1);
t192 = dph3.*t147;
t182 = t179+t180+t181-t192;
t183 = t150+t151+t152-t160;
t184 = t155+t157+t158-t161;
t185 = t6.*t126.*t127.*(1.0./2.0);
t186 = cos(th4);
t187 = cos(ph4);
t188 = sin(ph4);
t189 = sin(th4);
t193 = dth4.^2;
t194 = dph4.*t187.*t189.*(3.0./1.0e1);
t195 = dth4.*t186.*t188.*(3.0./1.0e1);
t196 = t194+t195;
t197 = dph4.*t186.*t188.*(3.0./1.0e1);
t198 = dth4.*t187.*t189.*(3.0./1.0e1);
t199 = t197+t198;
t200 = dph4.*t188.*t189.*(3.0./1.0e1);
t210 = dth4.*t186.*t187.*(3.0./1.0e1);
t201 = t200-t210;
t202 = dph4.*t186.*t187.*(3.0./1.0e1);
t208 = dth4.*t188.*t189.*(3.0./1.0e1);
t203 = t202-t208;
t204 = dph4.*t199;
t205 = dth4.*t196;
t206 = ddth4.*t188.*t189.*(3.0./1.0e1);
t214 = ddph4.*t186.*t187.*(3.0./1.0e1);
t207 = t41+t42+t43+t60+t61+t62-t102-t104+t150+t151+t152-t160+t204+t205+t206-t214;
t209 = dph4.*t203;
t211 = ddph4.*t186.*t188.*(3.0./1.0e1);
t212 = ddth4.*t187.*t189.*(3.0./1.0e1);
t215 = dth4.*t201;
t213 = t37+t39+t40+t55+t57+t58-t63-t103+t155+t157+t158-t161+t209+t211+t212-t215;
t216 = dph4.*t188.*t189.*(3.0./5.0);
t217 = dph4.*t187.*t189.*(3.0./5.0);
t218 = dth4.*t186.*t188.*(3.0./5.0);
t219 = t217+t218;
t220 = dph4.*t186.*t187.*(3.0./5.0);
t221 = t220-dth4.*t188.*t189.*(3.0./5.0);
t222 = dph4.*t186.*t188.*(3.0./5.0);
t223 = dth4.*t187.*t189.*(3.0./5.0);
t224 = t222+t223;
t225 = ddth4.*t189.*(3.0./1.0e1);
t226 = t186.*t193.*(3.0./1.0e1);
t227 = t225+t226;
t228 = ddth4.*t186.*(3.0./1.0e1);
t229 = dph4.*t196;
t230 = dth4.*t199;
t231 = ddph4.*t188.*t189.*(3.0./1.0e1);
t232 = t229+t230+t231-ddth4.*t186.*t187.*(3.0./1.0e1);
t233 = dth4.*t203;
t234 = ddph4.*t187.*t189.*(3.0./1.0e1);
t235 = ddth4.*t186.*t188.*(3.0./1.0e1);
t236 = t233+t234+t235-dph4.*t201;
t237 = t204+t205+t206-t214;
t238 = t209+t211+t212-t215;
t239 = t6.*t186.*t187.*(1.0./2.0);
H = reshape([0.0,0.0,t2.*t3,t2.*t5,-t4,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t5,t3,0.0,0.0,0.0,0.0,-dth1.*t4.*t23-t20.*t22.*t67+t20.*t21.*(t68-dth1.*t2.*t3),t67.*t75+t80.*(t68-dth1.*t2.*t3)+dth1.*t4.*t20.*t70,0.0,0.0,0.0,0.0,0.0,0.0,-dth1.*t4.*t130-t67.*t126.*t129+t126.*t127.*(t68-dth1.*t2.*t3),dth1.*t4.*t126-t67.*t129.*t130+t127.*t130.*(t68-dth1.*t2.*t3),0.0,0.0,0.0,0.0,0.0,-dth1.*t4.*t189-t67.*t186.*t188+t186.*t187.*(t68-dth1.*t2.*t3),dth1.*t4.*t186-t67.*t188.*t189+t187.*t189.*(t68-dth1.*t2.*t3),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,-t20.*t22.*t73-t20.*t21.*t77,t73.*t75-t77.*t80,0.0,0.0,0.0,0.0,0.0,0.0,-t73.*t126.*t129-t77.*t126.*t127,-t73.*t129.*t130-t77.*t127.*t130,0.0,0.0,0.0,0.0,0.0,-t73.*t186.*t188-t77.*t186.*t187,-t73.*t188.*t189-t77.*t187.*t189,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t20.*t21,t20.*t22,-t23,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,dth2.*t23.^2.*(-3.0./5.0)-t20.*t22.*t84+t20.*t21.*(t81-dth2.*t20.*t21.*(3.0./5.0)),t75.*t84+t80.*(t81-dth2.*t20.*t21.*(3.0./5.0))+dth2.*t20.*t23.*t70.*(3.0./5.0),-t22,t21,0.0,0.0,0.0,0.0,dth2.*t23.*t130.*(-3.0./5.0)-t84.*t126.*t129+t126.*t127.*(t81-dth2.*t20.*t21.*(3.0./5.0)),dth2.*t23.*t126.*(3.0./5.0)-t84.*t129.*t130+t127.*t130.*(t81-dth2.*t20.*t21.*(3.0./5.0)),0.0,0.0,0.0,0.0,0.0,dth2.*t23.*t189.*(-3.0./5.0)-t84.*t186.*t188+t186.*t187.*(t81-dth2.*t20.*t21.*(3.0./5.0)),dth2.*t23.*t186.*(3.0./5.0)-t84.*t188.*t189+t187.*t189.*(t81-dth2.*t20.*t21.*(3.0./5.0)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t20.*t21.*t87-t20.*t22.*t90,t75.*t90-t80.*t87,0.0,0.0,1.0,0.0,0.0,0.0,-t87.*t126.*t127-t90.*t126.*t129,-t87.*t127.*t130-t90.*t129.*t130,0.0,0.0,0.0,0.0,0.0,-t87.*t186.*t187-t90.*t186.*t188,-t87.*t187.*t189-t90.*t188.*t189,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,dth3.*t130.^2.*(-1.8e1./2.5e1)-t126.*t129.*t165+t126.*t127.*(t162-dth3.*t126.*t127.*(1.8e1./2.5e1)),dth3.*t126.*t130.*(1.8e1./2.5e1)-t129.*t130.*t165+t127.*t130.*(t162-dth3.*t126.*t127.*(1.8e1./2.5e1)),t127,0.0,0.0,0.0,0.0,dth3.*t130.*t189.*(-1.8e1./2.5e1)-t165.*t186.*t188+t186.*t187.*(t162-dth3.*t126.*t127.*(1.8e1./2.5e1)),dth3.*t130.*t186.*(1.8e1./2.5e1)-t165.*t188.*t189+t187.*t189.*(t162-dth3.*t126.*t127.*(1.8e1./2.5e1)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t126.*t127.*t167-t126.*t129.*t170,-t127.*t130.*t167-t129.*t130.*t170,0.0,1.0,0.0,0.0,0.0,-t167.*t186.*t187-t170.*t186.*t188,-t167.*t187.*t189-t170.*t188.*t189,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,dth4.*t189.^2.*(-3.0./5.0)-t186.*t188.*t219+t186.*t187.*(t216-dth4.*t186.*t187.*(3.0./5.0)),dth4.*t186.*t189.*(3.0./5.0)-t188.*t189.*t219+t187.*t189.*(t216-dth4.*t186.*t187.*(3.0./5.0)),t187,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t186.*t187.*t221-t186.*t188.*t224,-t187.*t189.*t221-t188.*t189.*t224,0.0,1.0,0.0,0.0,0.0,0.0,t2.*t7.*(-4.9e1./5.0),0.0,0.0,0.0,0.0,t17-t6.*t19.*(1.0./2.0),t6.*t11.*(-1.0./2.0)+t6.*t13.*(1.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2.*(-4.9e1./5.0),t4.*t8.*(-4.9e1./5.0),-dsi1.*t3.*t4,-dsi1.*t4.*t5,-dsi1.*t2,t3.*t4.*t6.*(-1.0./2.0)-t4.*t5.*t6.*(1.0./2.0),t2.*t3.*t6.*t7.*(1.0./2.0)+t2.*t5.*t6.*t7.*(1.0./2.0),t2.*t3.*t6.*t8.*(1.0./2.0)+t2.*t5.*t6.*t8.*(1.0./2.0),-t23.*t101+t20.*t21.*t94-t20.*t22.*t98,t75.*t98+t80.*t94+t20.*t70.*t101,0.0,0.0,0.0,0.0,0.0,0.0,-t101.*t130+t94.*t126.*t127-t98.*t126.*t129,t101.*t126+t94.*t127.*t130-t98.*t129.*t130,0.0,0.0,0.0,0.0,0.0,-t101.*t189+t94.*t186.*t187-t98.*t186.*t188,t101.*t186+t94.*t187.*t189-t98.*t188.*t189,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-dth1.*t3-dsi1.*t2.*t5,-dth1.*t5+dsi1.*t2.*t3,0.0,t2.*t3.*t6.*(1.0./2.0)-t2.*t5.*t6.*(1.0./2.0),t6.*t11.*(-1.0./2.0)-t6.*t13.*(1.0./2.0),t17+t6.*t19.*(1.0./2.0),-t20.*t22.*(t37+t39+t40-dth1.*t33)+t20.*t21.*(t41+t42+t43-ddph1.*t2.*t3.*(1.0./2.0)),t75.*t135+t80.*t134,0.0,0.0,0.0,0.0,0.0,0.0,t126.*t127.*t134-t126.*t129.*t135,t127.*t130.*t134-t129.*t130.*t135,0.0,0.0,0.0,0.0,0.0,t134.*t186.*t187-t135.*t186.*t188,t134.*t187.*t189-t135.*t188.*t189,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t20.*t69.*(-4.9e1./5.0)-t59.*t123+t64.*t121+t20.*t69.*(t105+t106-t4.*t34.*(1.0./2.0)-t23.*t35.*(3.0./1.0e1)),0.0,0.0,0.0,0.0,t125-t6.*t75.*(1.0./2.0),t6.*t121.*(-1.0./2.0)+t6.*t123.*(1.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t20.*(-4.9e1./5.0)-t23.*t117+t20.*(t105+t106-t4.*t34.*(1.0./2.0)-t23.*t35.*(3.0./1.0e1))+t21.*t23.*t59+t22.*t23.*t64+t20.*t21.*t110-t20.*t22.*t114,t23.*t70.*(-4.9e1./5.0)+t75.*t114+t80.*t110+t23.*t70.*(t105+t106-t4.*t34.*(1.0./2.0)-t23.*t35.*(3.0./1.0e1))+t20.*t70.*t117-t20.*t21.*t59.*t70-t20.*t22.*t64.*t70,-dsi2.*t21.*t23,-dsi2.*t22.*t23,-dsi2.*t20,t6.*t21.*t23.*(-1.0./2.0)-t6.*t22.*t23.*(1.0./2.0),t6.*t20.*t21.*t69.*(1.0./2.0)+t6.*t20.*t22.*t69.*(1.0./2.0),t6.*t20.*t21.*t70.*(1.0./2.0)+t6.*t20.*t22.*t70.*(1.0./2.0),-t117.*t130+t110.*t126.*t127-t114.*t126.*t129,t117.*t126+t110.*t127.*t130-t114.*t129.*t130,0.0,0.0,0.0,0.0,0.0,-t117.*t189+t110.*t186.*t187-t114.*t186.*t188,t117.*t186+t110.*t187.*t189-t114.*t188.*t189,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t20.*t21.*(t60+t61+t62-ddph2.*t20.*t21.*(3.0./1.0e1))+t20.*t22.*t59-t20.*t21.*t64-t20.*t22.*t118,-t59.*t75-t64.*t80+t75.*t118+t80.*t138,-dth2.*t21-dsi2.*t20.*t22,-dth2.*t22+dsi2.*t20.*t21,0.0,t6.*t20.*t21.*(1.0./2.0)-t6.*t20.*t22.*(1.0./2.0),t6.*t121.*(-1.0./2.0)-t6.*t123.*(1.0./2.0),t125+t6.*t75.*(1.0./2.0),-t118.*t126.*t129+t126.*t127.*t138,-t118.*t129.*t130+t127.*t130.*t138,0.0,0.0,0.0,0.0,0.0,-t118.*t186.*t188+t138.*t186.*t187,-t118.*t188.*t189+t138.*t187.*t189,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t126.*(-4.9e1./5.0)-t130.*t173+t126.*(t105+t106+t174-t4.*t34.*(1.0./2.0)-t23.*t35.*(3.0./1.0e1)-t130.*t139.*(9.0./2.5e1))+t129.*t130.*t153+t127.*t130.*t159+t126.*t127.*t178-t126.*t129.*t182,t130.*(-4.9e1./5.0)+t126.*t173+t130.*(t105+t106+t174-t4.*t34.*(1.0./2.0)-t23.*t35.*(3.0./1.0e1)-t130.*t139.*(9.0./2.5e1))-t126.*t129.*t153-t126.*t127.*t159+t127.*t130.*t178-t129.*t130.*t182,0.0,0.0,t6.*t127.*t130.*(-1.0./2.0)-t6.*t129.*t130.*(1.0./2.0),0.0,t185+t6.*t126.*t129.*(1.0./2.0),-t173.*t189+t178.*t186.*t187-t182.*t186.*t188,t173.*t186+t178.*t187.*t189-t182.*t188.*t189,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t126.*t127.*t153+t126.*t129.*t159+t126.*t127.*t183-t126.*t129.*t184,-t127.*t130.*t153+t129.*t130.*t159+t127.*t130.*t183-t129.*t130.*t184,-dth3.*t129,0.0,t185-t6.*t126.*t129.*(1.0./2.0),t6.*t127.*(-1.0./2.0)-t6.*t129.*(1.0./2.0),t6.*t127.*t130.*(1.0./2.0)-t6.*t129.*t130.*(1.0./2.0),t183.*t186.*t187-t184.*t186.*t188,t183.*t187.*t189-t184.*t188.*t189,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t186.*(-4.9e1./5.0)+t186.*(t105+t106+t174+t228-t4.*t34.*(1.0./2.0)-t23.*t35.*(3.0./1.0e1)-t130.*t139.*(9.0./2.5e1)-t189.*t193.*(3.0./1.0e1))-t189.*t227+t188.*t189.*t207+t187.*t189.*t213+t186.*t187.*t232-t186.*t188.*t236,t189.*(-4.9e1./5.0)+t189.*(t105+t106+t174+t228-t4.*t34.*(1.0./2.0)-t23.*t35.*(3.0./1.0e1)-t130.*t139.*(9.0./2.5e1)-t189.*t193.*(3.0./1.0e1))+t186.*t227-t186.*t188.*t207-t186.*t187.*t213+t187.*t189.*t232-t188.*t189.*t236,0.0,0.0,t6.*t187.*t189.*(-1.0./2.0)-t6.*t188.*t189.*(1.0./2.0),0.0,t239+t6.*t186.*t188.*(1.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t186.*t187.*t207+t186.*t188.*t213+t186.*t187.*t237-t186.*t188.*t238,-t187.*t189.*t207+t188.*t189.*t213+t187.*t189.*t237-t188.*t189.*t238,-dth4.*t188,0.0,t239-t6.*t186.*t188.*(1.0./2.0),t6.*t187.*(-1.0./2.0)-t6.*t188.*(1.0./2.0),t6.*t187.*t189.*(1.0./2.0)-t6.*t188.*t189.*(1.0./2.0)],[30, 20]);
