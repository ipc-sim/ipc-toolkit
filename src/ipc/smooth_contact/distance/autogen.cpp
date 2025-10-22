#include "autogen.hpp"

#include <cmath>

namespace ipc::autogen {

// hess is (144×1) flattened in column-major order
void edge_edge_closest_point_hessian_a(
    double ea0_x,
    double ea0_y,
    double ea0_z,
    double ea1_x,
    double ea1_y,
    double ea1_z,
    double eb0_x,
    double eb0_y,
    double eb0_z,
    double eb1_x,
    double eb1_y,
    double eb1_z,
    double hess[144])
{
    const auto t0 = ea0_y * eb0_y;
    const auto t1 = ea0_x * eb0_x;
    const auto t2 = 2 * t1;
    const auto t3 = ea0_y * eb1_y;
    const auto t4 = ea0_x * eb1_x;
    const auto t5 = 2 * t4;
    const auto t6 = ea0_z * eb0_z;
    const auto t7 = ea0_z * eb1_z;
    const auto t8 = ea0_x * eb0_y;
    const auto t9 = ea1_x * eb1_y;
    const auto t10 = ea0_x * eb0_z;
    const auto t11 = ea1_x * eb1_z;
    const auto t12 = ea1_y * eb0_y;
    const auto t13 = ea1_y * eb1_y;
    const auto t14 = ea1_z * eb0_z;
    const auto t15 = ea1_z * eb1_z;
    const auto t16 = 2 * t0;
    const auto t17 = 2 * t3;
    const auto t18 = ea1_x * eb0_x;
    const auto t19 = ea1_x * eb1_x;
    const auto t20 = ea0_y * eb1_x;
    const auto t21 = ea1_y * eb0_x;
    const auto t22 = ea0_y * eb0_z;
    const auto t23 = ea1_y * eb1_z;
    const auto t24 = 2 * t6;
    const auto t25 = 2 * t7;
    const auto t26 = ea0_z * eb1_x;
    const auto t27 = ea1_z * eb0_x;
    const auto t28 = ea0_z * eb1_y;
    const auto t29 = ea1_z * eb0_y;
    const auto t30 = 2 * t18;
    const auto t31 = 2 * t19;
    const auto t32 = 2 * t12;
    const auto t33 = 2 * t13;
    const auto t34 = std::pow(ea0_x, 2);
    const auto t35 = std::pow(eb0_y, 2);
    const auto t36 = std::pow(eb0_z, 2);
    const auto t37 = std::pow(eb1_y, 2);
    const auto t38 = std::pow(eb1_z, 2);
    const auto t39 = std::pow(ea0_y, 2);
    const auto t40 = std::pow(eb0_x, 2);
    const auto t41 = std::pow(eb1_x, 2);
    const auto t42 = std::pow(ea0_z, 2);
    const auto t43 = std::pow(ea1_x, 2);
    const auto t44 = std::pow(ea1_y, 2);
    const auto t45 = std::pow(ea1_z, 2);
    const auto t46 = ea0_x * t35;
    const auto t47 = 2 * ea1_x;
    const auto t48 = ea0_x * t36;
    const auto t49 = ea0_x * t37;
    const auto t50 = ea0_x * t38;
    const auto t51 = 2 * eb0_y;
    const auto t52 = eb1_y * t34;
    const auto t53 = 2 * eb0_z;
    const auto t54 = eb1_z * t34;
    const auto t55 = 2 * ea0_y;
    const auto t56 = ea1_y * t40;
    const auto t57 = ea1_y * t36;
    const auto t58 = ea1_y * t41;
    const auto t59 = ea1_y * t38;
    const auto t60 = 2 * eb0_x;
    const auto t61 = eb1_x * t39;
    const auto t62 = eb1_z * t39;
    const auto t63 = 2 * ea0_z;
    const auto t64 = ea1_z * t40;
    const auto t65 = ea1_z * t35;
    const auto t66 = ea1_z * t41;
    const auto t67 = ea1_z * t37;
    const auto t68 = eb1_x * t42;
    const auto t69 = eb1_y * t42;
    const auto t70 = eb1_y * t43;
    const auto t71 = eb1_z * t43;
    const auto t72 = eb1_x * t44;
    const auto t73 = eb1_z * t44;
    const auto t74 = eb1_x * t45;
    const auto t75 = eb1_y * t45;
    const auto t76 = -t0 * t2 + t0 * t5 + 4 * t10 * t11 + t12 * t2 + t12 * t24
        - t12 * t25 - t12 * t30 + t12 * t31 - t12 * t5 - t13 * t2 - t13 * t24
        + t13 * t25 + t13 * t30 - t13 * t31 + t13 * t5 + t14 * t16 - t14 * t17
        + t14 * t2 - t14 * t30 + t14 * t31 - t14 * t32 + t14 * t33 - t14 * t5
        - t15 * t16 + t15 * t17 - t15 * t2 + t15 * t30 - t15 * t31 + t15 * t32
        - t15 * t33 + t15 * t5 + t16 * t18 - t16 * t19 - t16 * t6 + t16 * t7
        - t17 * t18 + t17 * t19 + t17 * t6 - t17 * t7 + t18 * t24 - t18 * t25
        - t19 * t24 + t19 * t25 + t2 * t3 - t2 * t6 + t2 * t7 + 4 * t20 * t21
        + 4 * t22 * t23 + 4 * t26 * t27 + 4 * t28 * t29 - t3 * t5 + t34 * t35
        + t34 * t36 + t34 * t37 + t34 * t38 + t35 * t42 + t35 * t43 + t35 * t45
        + t36 * t39 + t36 * t43 + t36 * t44 + t37 * t42 + t37 * t43 + t37 * t45
        + t38 * t39 + t38 * t43 + t38 * t44 + t39 * t40 + t39 * t41 + t40 * t42
        + t40 * t44 + t40 * t45 + t41 * t42 + t41 * t44 + t41 * t45 - t46 * t47
        - t47 * t48 - t47 * t49 - t47 * t50 + t5 * t6 - t5 * t7 - t51 * t52
        - t51 * t69 - t51 * t70 - t51 * t75 - t53 * t54 - t53 * t62 - t53 * t71
        - t53 * t73 - t55 * t56 - t55 * t57 - t55 * t58 - t55 * t59 - t60 * t61
        - t60 * t68 - t60 * t72 - t60 * t74 - t63 * t64 - t63 * t65 - t63 * t66
        - t63 * t67 + 4 * t8 * t9;
    const auto t77 = 1.0 / t76;
    const auto t78 = eb1_z * t53;
    const auto t79 = t36 + t38 - t78;
    const auto t80 = eb1_y * t51;
    const auto t81 = t35 + t37 - t80;
    const auto t82 = t79 + t81;
    const auto t83 = ea0_x - eb0_x;
    const auto t84 = eb0_x - eb1_x;
    const auto t85 = t83 * t84;
    const auto t86 = ea0_y - eb0_y;
    const auto t87 = eb0_y - eb1_y;
    const auto t88 = t86 * t87;
    const auto t89 = ea0_z - eb0_z;
    const auto t90 = eb0_z - eb1_z;
    const auto t91 = t89 * t90;
    const auto t92 = t88 + t91;
    const auto t93 = t85 + t92;
    const auto t94 = -t14 + t15 + t6 - t7;
    const auto t95 = t0 - t12 + t13 - t3;
    const auto t96 = t94 + t95;
    const auto t97 = t1 - t18 + t19 - t4;
    const auto t98 = t96 + t97;
    const auto t99 = t93 * t98;
    const auto t100 = t82 * t99;
    const auto t101 = ea0_x - ea1_x;
    const auto t102 = t101 * t83;
    const auto t103 = ea0_y - ea1_y;
    const auto t104 = t103 * t86;
    const auto t105 = ea0_z - ea1_z;
    const auto t106 = t105 * t89;
    const auto t107 = t102 + t104 + t106;
    const auto t108 = eb1_x * t60;
    const auto t109 = -t108 + t40 + t41;
    const auto t110 = t109 + t82;
    const auto t111 = std::pow(t76, -2);
    const auto t112 = 2 * ea0_x;
    const auto t113 = eb1_y * t112;
    const auto t114 = eb1_z * t112;
    const auto t115 = -ea1_x * t35 - ea1_x * t36 - ea1_x * t37 - ea1_x * t38
        - eb0_x * t0 + eb0_x * t12 - eb0_x * t13 + eb0_x * t14 - eb0_x * t15
        + eb0_x * t3 - eb0_x * t6 + eb0_x * t7 - eb0_y * t113 - eb0_z * t114
        + eb1_x * t0 - eb1_x * t12 + eb1_x * t13 - eb1_x * t14 + eb1_x * t15
        - eb1_x * t3 + eb1_x * t6 - eb1_x * t7 + t11 * t53 + t46 + t48 + t49
        + t50 + t51 * t9;
    const auto t116 = std::pow(t115, 2);
    const auto t117 = t111 * t116;
    const auto t118 = 4 * t117;
    const auto t119 = 2 * t77;
    const auto t120 = t115 * t93;
    const auto t121 = t120 * t84;
    const auto t122 = t107 * t110;
    const auto t123 = t115 * t84;
    const auto t124 = t119 * t98;
    const auto t125 = t110 * t115;
    const auto t126 = ea1_x + eb0_x - t112;
    const auto t127 = t119 * t126;
    const auto t128 = t110 + t123 * t124 + t125 * t127 - std::pow(t84, 2);
    const auto t129 = t84 * t87;
    const auto t130 =
        eb0_x * eb0_y - eb0_x * eb1_y - eb0_y * eb1_x + eb1_x * eb1_y;
    const auto t131 = -t101 * t83 - t103 * t86 - t105 * t89;
    const auto t132 = t110 * t131;
    const auto t133 = t132 * t77;
    const auto t134 = t130 * t99;
    const auto t135 = eb0_x * t55;
    const auto t136 = eb1_z * t55;
    const auto t137 = ea1_y * eb1_x;
    const auto t138 = -ea0_y * t36 - ea0_y * t38 - ea0_y * t40 - ea0_y * t41
        + eb0_y * t1 - eb0_y * t14 + eb0_y * t15 - eb0_y * t18 + eb0_y * t19
        - eb0_y * t4 + eb0_y * t6 - eb0_y * t7 + eb0_z * t136 + eb1_x * t135
        - eb1_y * t1 + eb1_y * t14 - eb1_y * t15 + eb1_y * t18 - eb1_y * t19
        + eb1_y * t4 - eb1_y * t6 + eb1_y * t7 - t137 * t60 - t23 * t53 + t56
        + t57 + t58 + t59;
    const auto t139 = t110 * t138;
    const auto t140 = t126 * t77;
    const auto t141 = ea1_y + eb0_y - t55;
    const auto t142 = t84 * t93;
    const auto t143 = t138 * t142;
    const auto t144 = t77 * t98;
    const auto t145 = t144 * t84;
    const auto t146 = t119
        * (4 * t110 * t111 * t115 * t131 * t138 + t110 * t115 * t141 * t77
           + 4 * t111 * t115 * t138 * t93 * t98 + t115 * t77 * t87 * t93
           + t115 * t77 * t87 * t98 - t129 - t130 * t133 - t134 * t77
           - t138 * t145 - t139 * t140 - t143 * t77);
    const auto t147 = t84 * t90;
    const auto t148 =
        eb0_x * eb0_z - eb0_x * eb1_z - eb0_z * eb1_x + eb1_x * eb1_z;
    const auto t149 = t148 * t99;
    const auto t150 = eb0_x * t63;
    const auto t151 = eb0_y * t63;
    const auto t152 = ea1_z * eb1_x;
    const auto t153 = ea1_z * eb1_y;
    const auto t154 = -ea0_z * t35 - ea0_z * t37 - ea0_z * t40 - ea0_z * t41
        + eb0_z * t0 + eb0_z * t1 - eb0_z * t12 + eb0_z * t13 - eb0_z * t18
        + eb0_z * t19 - eb0_z * t3 - eb0_z * t4 + eb1_x * t150 + eb1_y * t151
        - eb1_z * t0 - eb1_z * t1 + eb1_z * t12 - eb1_z * t13 + eb1_z * t18
        - eb1_z * t19 + eb1_z * t3 + eb1_z * t4 - t152 * t60 - t153 * t51 + t64
        + t65 + t66 + t67;
    const auto t155 = t110 * t154;
    const auto t156 = ea1_z + eb0_z - t63;
    const auto t157 = t142 * t154;
    const auto t158 = t119
        * (4 * t110 * t111 * t115 * t131 * t154 + t110 * t115 * t156 * t77
           + 4 * t111 * t115 * t154 * t93 * t98 + t115 * t77 * t90 * t93
           + t115 * t77 * t90 * t98 - t133 * t148 - t140 * t155 - t145 * t154
           - t147 - t149 * t77 - t157 * t77);
    const auto t159 = 4 * t77;
    const auto t160 = t77
        * (-t100 * t119 + 2 * t107 * t110 * t77 * t82
           + 2 * t110 * t115 * t77 * t83 + 8 * t111 * t116 * t93 * t98
           - 8 * t117 * t122 - t121 * t159 - t128);
    const auto t161 = t122 * t130;
    const auto t162 = t125 * t86;
    const auto t163 = t120 * t87;
    const auto t164 = t119 * t163;
    const auto t165 = t124 * t84;
    const auto t166 = 8 * t111;
    const auto t167 = t138 * t166;
    const auto t168 = t115 * t122;
    const auto t169 = t115 * t99;
    const auto t170 = t167 * t169;
    const auto t171 = t77
        * (t119 * t134 + t119 * t143 - t119 * t161 + t119 * t162 + t127 * t139
           + t129 + t138 * t165 - t164 + t167 * t168 - t170);
    const auto t172 = t122 * t148;
    const auto t173 = t125 * t89;
    const auto t174 = t120 * t90;
    const auto t175 = t119 * t174;
    const auto t176 = t154 * t166;
    const auto t177 = t169 * t176;
    const auto t178 = t77
        * (t119 * t149 + t119 * t157 - t119 * t172 + t119 * t173 + t127 * t155
           + t147 + t154 * t165 + t168 * t176 - t175 - t177);
    const auto t179 = ea0_x + eb1_x - t60;
    const auto t180 = t131 * t159;
    const auto t181 = t101 * t119;
    const auto t182 = t115 * t124;
    const auto t183 = t119 * t96;
    const auto t184 = t125 * t131;
    const auto t185 = ea0_x * t0 - ea0_x * t12 + ea0_x * t13 - ea0_x * t14
        + ea0_x * t15 - ea0_x * t3 + ea0_x * t6 - ea0_x * t7 - ea1_x * t0
        + ea1_x * t12 - ea1_x * t13 + ea1_x * t14 - ea1_x * t15 + ea1_x * t3
        - ea1_x * t6 + ea1_x * t7 - eb0_x * t39 - eb0_x * t42 - eb0_x * t44
        - eb0_x * t45 - t137 * t55 - t152 * t63 + t21 * t55 + t27 * t63 + t61
        + t68 + t72 + t74;
    const auto t186 = 8 * t185;
    const auto t187 = t111 * t186;
    const auto t188 = t110 - t123 * t180 - t125 * t181 + t132 * t183
        - t179 * t182 + t179 * t84 - t184 * t187;
    const auto t189 = t101 * t84;
    const auto t190 = 2 * t126;
    const auto t191 = t110 * t127;
    const auto t192 = t119 * t142;
    const auto t193 = -t120 * t181 - t169 * t187 + t183 * t99 + t185 * t192;
    const auto t194 =
        t165 * t185 + t185 * t191 + t189 + t190 * t84 + t193 + t98;
    const auto t195 = t77 * (-t188 - t194 - t93);
    const auto t196 = -ea0_y * eb0_x - ea1_x * t51 + eb0_y * t112 - t113 - t137
        + t20 + t21 + 2 * t9;
    const auto t197 = t119 * t196;
    const auto t198 = t115 * t180;
    const auto t199 = ea1_x * eb0_y;
    const auto t200 = -ea0_y * t1 + ea0_y * t14 - ea0_y * t15 + ea0_y * t18
        - ea0_y * t19 + ea0_y * t4 - ea0_y * t6 + ea0_y * t7 + ea1_y * t1
        - ea1_y * t14 + ea1_y * t15 - ea1_y * t18 + ea1_y * t19 - ea1_y * t4
        + ea1_y * t6 - ea1_y * t7 + eb0_y * t34 + eb0_y * t42 + eb0_y * t43
        + eb0_y * t45 - t112 * t199 + t112 * t9 + t153 * t63 - t29 * t63 - t52
        - t69 - t70 - t75;
    const auto t201 = t166 * t200;
    const auto t202 = t132 * t197 - t184 * t201 + t198 * t87;
    const auto t203 = t103 * t119;
    const auto t204 = t120 * t203;
    const auto t205 = t192 * t200;
    const auto t206 = t197 * t99;
    const auto t207 = -t169 * t201;
    const auto t208 = ea0_y + eb1_y - t51;
    const auto t209 =
        t125 * t203 + t182 * t208 + t204 + t205 + t206 + t207 - t208 * t84;
    const auto t210 = t103 * t84;
    const auto t211 = t165 * t200 - t190 * t87 + t191 * t200 - t210;
    const auto t212 = t77 * (t202 + t209 + t211);
    const auto t213 = -ea0_z * eb0_x - ea1_x * t53 + eb0_z * t112 + 2 * t11
        - t114 - t152 + t26 + t27;
    const auto t214 = t119 * t213;
    const auto t215 = ea1_x * eb0_z;
    const auto t216 = ea1_y * eb0_z;
    const auto t217 = -ea0_z * t0 - ea0_z * t1 + ea0_z * t12 - ea0_z * t13
        + ea0_z * t18 - ea0_z * t19 + ea0_z * t3 + ea0_z * t4 + ea1_z * t0
        + ea1_z * t1 - ea1_z * t12 + ea1_z * t13 - ea1_z * t18 + ea1_z * t19
        - ea1_z * t3 - ea1_z * t4 + eb0_z * t34 + eb0_z * t39 + eb0_z * t43
        + eb0_z * t44 + t11 * t112 - t112 * t215 - t216 * t55 + t23 * t55 - t54
        - t62 - t71 - t73;
    const auto t218 = t166 * t217;
    const auto t219 = t132 * t214 - t184 * t218 + t198 * t90;
    const auto t220 = t105 * t119;
    const auto t221 = t120 * t220;
    const auto t222 = t192 * t217;
    const auto t223 = t214 * t99;
    const auto t224 = -t169 * t218;
    const auto t225 = ea0_z + eb1_z - t53;
    const auto t226 =
        t125 * t220 + t182 * t225 + t221 + t222 + t223 + t224 - t225 * t84;
    const auto t227 = t105 * t84;
    const auto t228 = t165 * t217 - t190 * t90 + t191 * t217 - t227;
    const auto t229 = t77 * (t219 + t226 + t228);
    const auto t230 = t107 * t159;
    const auto t231 = -t122 * t183 + t123 * t230 + t168 * t187 - t182 * t83;
    const auto t232 = t77 * (t194 + t231 + 2 * t85 + t92);
    const auto t233 = t84 * t86;
    const auto t234 = t182 * t86 + t204 + t205 + t206 + t207 - t233;
    const auto t235 = t77
        * (2 * t107 * t110 * t196 * t77 + 4 * t107 * t115 * t77 * t87
           - t168 * t201 - t211 - t234);
    const auto t236 = t84 * t89;
    const auto t237 = t182 * t89 + t221 + t222 + t223 + t224 - t236;
    const auto t238 = t77
        * (2 * t107 * t110 * t213 * t77 + 4 * t107 * t115 * t77 * t90
           - t168 * t218 - t228 - t237);
    const auto t239 = std::pow(t87, 2);
    const auto t240 = t109 + t79;
    const auto t241 = t138 * t87;
    const auto t242 = t124 * t241;
    const auto t243 = t119 * t141;
    const auto t244 = t139 * t243;
    const auto t245 = std::pow(t138, 2);
    const auto t246 = 4 * t111 * t245;
    const auto t247 = t108 - t35 - t36 - t37 - t38 - t40 - t41 + t78 + t80;
    const auto t248 = t87 * t90;
    const auto t249 =
        eb0_y * eb0_z - eb0_y * eb1_z - eb0_z * eb1_y + eb1_y * eb1_z;
    const auto t250 = t249 * t99;
    const auto t251 = t141 * t155;
    const auto t252 = t139 * t156;
    const auto t253 = t154 * t87;
    const auto t254 = t253 * t93;
    const auto t255 = t138 * t90;
    const auto t256 = t255 * t93;
    const auto t257 = 4 * t111 * t138;
    const auto t258 = t131 * t155;
    const auto t259 = t154 * t99;
    const auto t260 = -t119
        * (t133 * t249 + t144 * t253 + t144 * t255 + t248 + t250 * t77
           + t251 * t77 + t252 * t77 + t254 * t77 + t256 * t77 + t257 * t258
           + t257 * t259);
    const auto t261 = t110 * t83;
    const auto t262 = t138 * t261;
    const auto t263 = t77
        * (2 * t110 * t130 * t131 * t77 - t119 * t262 - t125 * t243 + t129
           + 2 * t130 * t77 * t93 * t98 + 2 * t138 * t77 * t84 * t93 - t164
           - t167 * t184 - t170 - t182 * t87);
    const auto t264 = t139 * t86;
    const auto t265 = t119 * t132;
    const auto t266 = t240 * t99;
    const auto t267 = t77
        * (8 * t110 * t111 * t131 * t245 - t110 + 8 * t111 * t245 * t93 * t98
           - t119 * t264 - t119 * t266 + 4 * t138 * t77 * t87 * t93 + t239
           - t240 * t265 + t242 + t244);
    const auto t268 = t139 * t89;
    const auto t269 = t119 * t250 + t119 * t254 + t119 * t256 + t167 * t258
        + t167 * t259 + t248 + t249 * t265;
    const auto t270 = t77 * (t119 * t251 - t119 * t268 + t124 * t253 + t269);
    const auto t271 = t124 * t138;
    const auto t272 = -ea0_x * eb1_y + ea1_y * t60 + eb1_x * t55 - t135
        - 2 * t137 - t199 + t8 + t9;
    const auto t273 = t138 * t84;
    const auto t274 = t131 * t139;
    const auto t275 = t139 * t181 + t179 * t271 + t179 * t87 + t180 * t273
        + t187 * t274 + t265 * t272;
    const auto t276 = t101 * t87;
    const auto t277 = 2 * t141;
    const auto t278 = t185 * t87;
    const auto t279 = t110 * t243;
    const auto t280 = t138 * t93;
    const auto t281 = t119 * t93;
    const auto t282 = t119 * t99;
    const auto t283 = t138 * t187;
    const auto t284 = t181 * t280 + t272 * t282 + t278 * t281 + t283 * t99;
    const auto t285 = t124 * t278 + t185 * t279 + t276 + t277 * t84 + t284;
    const auto t286 = -t77 * (t275 + t285);
    const auto t287 = t94 + t97;
    const auto t288 =
        -8 * t110 * t111 * t131 * t138 * t200 + t180 * t241 + t265 * t287;
    const auto t289 = t110 + t139 * t203 + t208 * t271 + t208 * t87;
    const auto t290 = t200 * t87;
    const auto t291 = -t281 * t290;
    const auto t292 = t203 * t280;
    const auto t293 = t282 * t287;
    const auto t294 = t167 * t99;
    const auto t295 = -t200 * t294;
    const auto t296 = t103 * t87;
    const auto t297 = -t124 * t290 - t200 * t279 + t277 * t87 + t291 + t292
        + t293 + t295 + t296;
    const auto t298 = t93 + t98;
    const auto t299 = t77 * (-t288 - t289 - t297 - t298);
    const auto t300 = t217 * t87;
    const auto t301 = t281 * t300;
    const auto t302 = -t301;
    const auto t303 = -ea0_z * eb0_y - ea1_y * t53 + eb0_z * t55 - t136 - t153
        + 2 * t23 + t28 + t29;
    const auto t304 = t282 * t303;
    const auto t305 = -t304;
    const auto t306 = t220 * t280;
    const auto t307 = t217 * t294;
    const auto t308 = -t307;
    const auto t309 = t180 * t255;
    const auto t310 = t139 * t220 + t225 * t271 + t225 * t87;
    const auto t311 = t105 * t87;
    const auto t312 = -t124 * t300 - t217 * t279 + t277 * t90 + t311;
    const auto t313 = t77
        * (8 * t110 * t111 * t131 * t138 * t217 + 2 * t110 * t131 * t303 * t77
           - t302 - t305 - t306 - t308 - t309 - t310 - t312);
    const auto t314 = t83 * t87;
    const auto t315 = t119 * t122;
    const auto t316 =
        -t122 * t283 - t230 * t273 + t271 * t83 - t272 * t315 + t314;
    const auto t317 = t77 * (t285 + t316);
    const auto t318 = t271 * t86;
    const auto t319 = t122 * t167;
    const auto t320 = t200 * t319 - t230 * t241 - t287 * t315;
    const auto t321 = t85 + t98;
    const auto t322 = t77 * (t297 + t318 + t320 + t321 + 2 * t88 + t91);
    const auto t323 = t87 * t89;
    const auto t324 = t271 * t89;
    const auto t325 =
        t217 * t319 - t230 * t255 + t302 + t303 * t315 + t305 + t306 + t308;
    const auto t326 = t77 * (t312 + t323 + t324 + t325);
    const auto t327 = std::pow(t90, 2);
    const auto t328 = t109 + t81;
    const auto t329 = t154 * t90;
    const auto t330 = t124 * t329;
    const auto t331 = t119 * t156;
    const auto t332 = t155 * t331;
    const auto t333 = std::pow(t154, 2);
    const auto t334 = 4 * t111 * t333;
    const auto t335 = t154 * t261;
    const auto t336 = t77
        * (2 * t110 * t131 * t148 * t77 - t119 * t335 - t125 * t331 + t147
           + 2 * t148 * t77 * t93 * t98 + 2 * t154 * t77 * t84 * t93 - t175
           - t176 * t184 - t177 - t182 * t90);
    const auto t337 = t155 * t86;
    const auto t338 = t77 * (t119 * t252 - t119 * t337 + t124 * t255 + t269);
    const auto t339 = t155 * t89;
    const auto t340 = t328 * t99;
    const auto t341 = t77
        * (8 * t110 * t111 * t131 * t333 - t110 + 8 * t111 * t333 * t93 * t98
           - t119 * t339 - t119 * t340 + 4 * t154 * t77 * t90 * t93
           - t265 * t328 + t327 + t330 + t332);
    const auto t342 = t124 * t154;
    const auto t343 = -ea0_x * eb1_z + ea1_z * t60 + eb1_x * t63 + t10 + t11
        - t150 - 2 * t152 - t215;
    const auto t344 = t154 * t84;
    const auto t345 = t155 * t181 + t179 * t342 + t179 * t90 + t180 * t344
        + t187 * t258 + t265 * t343;
    const auto t346 = t101 * t90;
    const auto t347 = 2 * t156;
    const auto t348 = t185 * t90;
    const auto t349 = t110 * t331;
    const auto t350 = t154 * t93;
    const auto t351 = t181 * t350 + t187 * t259 + t281 * t348 + t282 * t343;
    const auto t352 = t124 * t348 + t185 * t349 + t346 + t347 * t84 + t351;
    const auto t353 = -t77 * (t345 + t352);
    const auto t354 = -ea0_y * eb1_z + ea1_z * t51 + eb1_y * t63 - t151
        - 2 * t153 - t216 + t22 + t23;
    const auto t355 =
        -8 * t110 * t111 * t131 * t154 * t200 + t180 * t253 + t265 * t354;
    const auto t356 = t155 * t203 + t208 * t342 + t208 * t90;
    const auto t357 = t200 * t90;
    const auto t358 = -t281 * t357;
    const auto t359 = t203 * t350;
    const auto t360 = t282 * t354;
    const auto t361 = t176 * t99;
    const auto t362 = -t200 * t361;
    const auto t363 = t103 * t90;
    const auto t364 = -t124 * t357 - t200 * t349 + t347 * t87 + t358 + t359
        + t360 + t362 + t363;
    const auto t365 = t77 * (-t355 - t356 - t364);
    const auto t366 = t95 + t97;
    const auto t367 =
        -8 * t110 * t111 * t131 * t154 * t217 + t180 * t329 + t265 * t366;
    const auto t368 = t110 + t155 * t220 + t225 * t342 + t225 * t90;
    const auto t369 = t217 * t90;
    const auto t370 = -t281 * t369;
    const auto t371 = t220 * t350;
    const auto t372 = t282 * t366;
    const auto t373 = -t217 * t361;
    const auto t374 = t105 * t90;
    const auto t375 = -t124 * t369 - t217 * t349 + t347 * t90 + t370 + t371
        + t372 + t373 + t374;
    const auto t376 = t77 * (-t298 - t367 - t368 - t375);
    const auto t377 = t83 * t90;
    const auto t378 = t122 * t187;
    const auto t379 =
        -t154 * t378 - t230 * t344 - t315 * t343 + t342 * t83 + t377;
    const auto t380 = t77 * (t352 + t379);
    const auto t381 = t86 * t90;
    const auto t382 = t342 * t86 + t381;
    const auto t383 = t122 * t176;
    const auto t384 = t200 * t383 - t230 * t253 - t315 * t354;
    const auto t385 = t77 * (t364 + t382 + t384);
    const auto t386 = t342 * t89;
    const auto t387 = t217 * t383 - t230 * t329 - t315 * t366;
    const auto t388 = t77 * (t321 + t375 + t386 + t387 + t88 + 2 * t91);
    const auto t389 = 4 * t116;
    const auto t390 = t77 * t99;
    const auto t391 = 2 * t111;
    const auto t392 = t122 * t159;
    const auto t393 = t115 * t138;
    const auto t394 = 4 * t390;
    const auto t395 = t391
        * (-t134 - t143 + t161 - t162 + t163 + t262 - t392 * t393
           + t393 * t394);
    const auto t396 = t115 * t154;
    const auto t397 = t391
        * (-t149 - t157 + t172 - t173 + t174 + t335 - t392 * t396
           + t394 * t396);
    const auto t398 = t119 * t261;
    const auto t399 = -t185 * t398 + t193 + t92;
    const auto t400 = t77 * (t188 + t399 - t85);
    const auto t401 = -t200 * t398 + t202 + 2 * t314;
    const auto t402 = t77 * (-t209 - t401);
    const auto t403 = -t217 * t398 + t219 + 2 * t377;
    const auto t404 = t77 * (-t226 - t403);
    const auto t405 = t77 * (-t231 - t399);
    const auto t406 = t77 * (t234 + t401);
    const auto t407 = t77 * (t237 + t403);
    const auto t408 = 2 * t93;
    const auto t409 = t138 * t154;
    const auto t410 = t391
        * (t122 * t249 - t250 - t254 - t256 + t268 + t337 + t392 * t409
           - t394 * t409);
    const auto t411 = t110 * t119;
    const auto t412 = t185 * t86;
    const auto t413 = -2 * t233 + t284 - t411 * t412;
    const auto t414 = t77 * (t275 + t413);
    const auto t415 = t411 * t86;
    const auto t416 = t200 * t415 + t291 + t292 + t293 + t295 + t85 + t91;
    const auto t417 = t77 * (t289 + t320 + t416 - t88);
    const auto t418 = 2 * t381;
    const auto t419 = t217 * t415;
    const auto t420 = t77 * (t310 + t325 - t418 + t419);
    const auto t421 = t77 * (-t316 - t413);
    const auto t422 = t77 * (-t288 - t318 - t416);
    const auto t423 = t77
        * (t218 * t274 + t265 * t303 + t301 + t304 - t306 + t307 - t309 - t323
           - t324 + t418 - t419);
    const auto t424 = t411 * t89;
    const auto t425 = -t185 * t424 - 2 * t236 + t351;
    const auto t426 = t77 * (t345 + t425);
    const auto t427 = t200 * t424 - 2 * t323 + t358 + t359 + t360 + t362;
    const auto t428 = t77 * (t356 + t384 + t427);
    const auto t429 = t217 * t424 + t370 + t371 + t372 + t373 + t85 + t88;
    const auto t430 = t77 * (t368 + t387 + t429 - t91);
    const auto t431 = t77 * (-t379 - t425);
    const auto t432 = t77 * (-t355 - t382 - t427);
    const auto t433 = t77 * (-t367 - t386 - t429);
    const auto t434 = -ea1_z * t63 + t42 + t45;
    const auto t435 = -ea1_y * t55 + t39 + t44;
    const auto t436 = t434 + t435;
    const auto t437 = t181 * t185;
    const auto t438 = t111 * std::pow(t185, 2);
    const auto t439 = 4 * t438;
    const auto t440 = -t110 * t131 * t436 * t77 + t131 + t132 * t439
        + t180 * t185 * t84 - t436 * t77 * t93 * t98 + t437 * t93 + t439 * t99;
    const auto t441 = t124 * t185;
    const auto t442 =
        -t0 - t1 + t12 - t13 + t14 - t15 + t18 - t19 + t3 + t4 - t6 + t7;
    const auto t443 = t101 * t179 + t110 * t437 + t179 * t441 + 2 * t189 + t442;
    const auto t444 =
        ea0_x * ea0_y - ea0_x * ea1_y - ea0_y * ea1_x + ea1_x * ea1_y;
    const auto t445 = -8 * t110 * t111 * t131 * t185 * t200
        - 4 * t131 * t200 * t77 * t84 + t180 * t278 + t265 * t444;
    const auto t446 = t181 * t200;
    const auto t447 = t124 * t179;
    const auto t448 = t103 * t179 - t110 * t446 - t200 * t447 + 2 * t276;
    const auto t449 = -t446 * t93;
    const auto t450 = t185 * t203;
    const auto t451 = t450 * t93;
    const auto t452 = t282 * t444;
    const auto t453 = -t187 * t200 * t99;
    const auto t454 = t101 * t208 + t110 * t450 + t208 * t441 + 2 * t210 + t449
        + t451 + t452 + t453;
    const auto t455 = t77 * (-t445 - t448 - t454);
    const auto t456 =
        ea0_x * ea0_z - ea0_x * ea1_z - ea0_z * ea1_x + ea1_x * ea1_z;
    const auto t457 = -8 * t110 * t111 * t131 * t185 * t217
        - 4 * t131 * t217 * t77 * t84 + t180 * t348 + t265 * t456;
    const auto t458 = t181 * t217;
    const auto t459 = t105 * t179 - t110 * t458 - t217 * t447 + 2 * t346;
    const auto t460 = -t458 * t93;
    const auto t461 = t185 * t220;
    const auto t462 = t461 * t93;
    const auto t463 = t282 * t456;
    const auto t464 = -t187 * t217 * t99;
    const auto t465 = t101 * t225 + t110 * t461 + t225 * t441 + 2 * t227 + t460
        + t462 + t463 + t464;
    const auto t466 = t77 * (-t457 - t459 - t465);
    const auto t467 = t441 * t83;
    const auto t468 = 2 * t104;
    const auto t469 = 2 * t106;
    const auto t470 = 8 * t438;
    const auto t471 = t159 * t93;
    const auto t472 = t77
        * (t101 * t185 * t471 - t102 - t107 * t186 * t77 * t84 - t122 * t470
           - t282 * t436 + t315 * t436 + t443 + t467 - t468 - t469
           + t470 * t99);
    const auto t473 = t101 * t86 + t124 * t412 + t449 + t451 + t452 + t453;
    const auto t474 = t230 * t84;
    const auto t475 = t200 * t378 + t200 * t474 - t230 * t278 - t315 * t444;
    const auto t476 = t77 * (t448 + t473 + t475);
    const auto t477 = t101 * t89 + t441 * t89 + t460 + t462 + t463 + t464;
    const auto t478 = t217 * t378 + t217 * t474 - t230 * t348 - t315 * t456;
    const auto t479 = t77 * (t459 + t477 + t478);
    const auto t480 = t103 * t208;
    const auto t481 = 2 * t296;
    const auto t482 = t110 * t200 * t203;
    const auto t483 = t124 * t208;
    const auto t484 = t200 * t483;
    const auto t485 = -ea1_x * t112 + t34 + t43;
    const auto t486 = t434 + t485;
    const auto t487 = std::pow(t200, 2);
    const auto t488 = t111 * t487;
    const auto t489 = 4 * t488;
    const auto t490 = -2 * t103 * t200 * t77 * t93 - t110 * t131 * t486 * t77
        - 4 * t131 * t200 * t77 * t87 + t131 + t132 * t489
        - t486 * t77 * t93 * t98 + t489 * t99;
    const auto t491 =
        ea0_y * ea0_z - ea0_y * ea1_z - ea0_z * ea1_y + ea1_y * ea1_z;
    const auto t492 = t201 * t217;
    const auto t493 = -4 * t131 * t200 * t77 * t90 - 4 * t131 * t217 * t77 * t87
        + t132 * t492 + t265 * t491;
    const auto t494 = t203 * t217;
    const auto t495 = t105 * t208 - t110 * t494 - t217 * t483 + 2 * t363;
    const auto t496 = -t494 * t93;
    const auto t497 = t200 * t220;
    const auto t498 = -t497 * t93;
    const auto t499 = t282 * t491;
    const auto t500 = t492 * t99;
    const auto t501 = t124 * t200;
    const auto t502 = t103 * t225 - t110 * t497 - t225 * t501 + 2 * t311 + t496
        + t498 + t499 + t500;
    const auto t503 = t77 * (-t493 - t495 - t502);
    const auto t504 = t103 * t83 - t501 * t83;
    const auto t505 = t77 * (t454 + t475 + t504);
    const auto t506 = t501 * t86;
    const auto t507 = 2 * t102 + t98;
    const auto t508 = t77
        * (-t103 * t200 * t471 - t104 + 2 * t107 * t110 * t486 * t77
           + 8 * t107 * t200 * t77 * t87 + 8 * t111 * t487 * t93 * t98
           - 8 * t122 * t488 - t282 * t486 - t469 + t480 + t481 - t482 - t484
           - t506 - t507);
    const auto t509 = t103 * t89 + t496 + t498 + t499 + t500 - t501 * t89;
    const auto t510 = -t122 * t492 + t230 * t300 + t230 * t357 - t315 * t491;
    const auto t511 = t77 * (t495 + t509 + t510);
    const auto t512 = t105 * t225;
    const auto t513 = 2 * t374;
    const auto t514 = t110 * t217 * t220;
    const auto t515 = t124 * t217;
    const auto t516 = t225 * t515;
    const auto t517 = t435 + t485;
    const auto t518 = std::pow(t217, 2);
    const auto t519 = t111 * t518;
    const auto t520 = 4 * t519;
    const auto t521 = -2 * t105 * t217 * t77 * t93 - t110 * t131 * t517 * t77
        - 4 * t131 * t217 * t77 * t90 + t131 + t132 * t520
        - t517 * t77 * t93 * t98 + t520 * t99;
    const auto t522 = t105 * t83 - t515 * t83;
    const auto t523 = t77 * (t465 + t478 + t522);
    const auto t524 = t105 * t86 - t515 * t86;
    const auto t525 = t77 * (t502 + t510 + t524);
    const auto t526 = t515 * t89;
    const auto t527 = t77
        * (-t105 * t217 * t471 - t106 + 2 * t107 * t110 * t517 * t77
           + 8 * t107 * t217 * t77 * t90 + 8 * t111 * t518 * t93 * t98
           - 8 * t122 * t519 - t282 * t517 - t468 - t507 + t512 + t513 - t514
           - t516 - t526);
    const auto t528 = t77 * (-t445 - t473 - t504);
    const auto t529 = t77 * (-t457 - t477 - t522);
    const auto t530 = t77 * (-t493 - t509 - t524);
    hess[0] = t119
        * (t100 * t77 - t107 * t110 * t77 * t82 + t118 * t122 - t118 * t99
           + t119 * t121 + t128);
    hess[1] = t146;
    hess[2] = t158;
    hess[3] = t160;
    hess[4] = t171;
    hess[5] = t178;
    hess[6] = t195;
    hess[7] = t212;
    hess[8] = t229;
    hess[9] = t232;
    hess[10] = t235;
    hess[11] = t238;
    hess[12] = t146;
    hess[13] = t119
        * (t110 * t131 * t240 * t77 - t119 * t241 * t93 - t132 * t246 - t239
           + t240 * t77 * t93 * t98 - t242 - t244 - t246 * t99 - t247);
    hess[14] = t260;
    hess[15] = t263;
    hess[16] = t267;
    hess[17] = t270;
    hess[18] = t286;
    hess[19] = t299;
    hess[20] = t313;
    hess[21] = t317;
    hess[22] = t322;
    hess[23] = t326;
    hess[24] = t158;
    hess[25] = t260;
    hess[26] = t119
        * (t110 * t131 * t328 * t77 - t132 * t334 - t247 - t281 * t329 - t327
           + t328 * t77 * t93 * t98 - t330 - t332 - t334 * t99);
    hess[27] = t336;
    hess[28] = t338;
    hess[29] = t341;
    hess[30] = t353;
    hess[31] = t365;
    hess[32] = t376;
    hess[33] = t380;
    hess[34] = t385;
    hess[35] = t388;
    hess[36] = t160;
    hess[37] = t263;
    hess[38] = t336;
    hess[39] = t391
        * (t110 * t131 * t82 - 2 * t115 * t261 + 2 * t115 * t84 * t93
           - t133 * t389 - t389 * t390 + t82 * t93 * t98);
    hess[40] = t395;
    hess[41] = t397;
    hess[42] = t400;
    hess[43] = t402;
    hess[44] = t404;
    hess[45] = t405;
    hess[46] = t406;
    hess[47] = t407;
    hess[48] = t171;
    hess[49] = t267;
    hess[50] = t338;
    hess[51] = t395;
    hess[52] = t391
        * (-t122 * t240 - t241 * t408 + t245 * t392 - t245 * t394 + 2 * t264
           + t266);
    hess[53] = t410;
    hess[54] = t414;
    hess[55] = t417;
    hess[56] = t420;
    hess[57] = t421;
    hess[58] = t422;
    hess[59] = t423;
    hess[60] = t178;
    hess[61] = t270;
    hess[62] = t341;
    hess[63] = t397;
    hess[64] = t410;
    hess[65] = t391
        * (-t122 * t328 - t329 * t408 + t333 * t392 - t333 * t394 + 2 * t339
           + t340);
    hess[66] = t426;
    hess[67] = t428;
    hess[68] = t430;
    hess[69] = t431;
    hess[70] = t432;
    hess[71] = t433;
    hess[72] = t195;
    hess[73] = t286;
    hess[74] = t353;
    hess[75] = t400;
    hess[76] = t414;
    hess[77] = t426;
    hess[78] = t119 * (-t440 - t443);
    hess[79] = t455;
    hess[80] = t466;
    hess[81] = t472;
    hess[82] = t476;
    hess[83] = t479;
    hess[84] = t212;
    hess[85] = t299;
    hess[86] = t365;
    hess[87] = t402;
    hess[88] = t417;
    hess[89] = t428;
    hess[90] = t455;
    hess[91] = t119 * (-t442 - t480 - t481 + t482 + t484 - t490);
    hess[92] = t503;
    hess[93] = t505;
    hess[94] = t508;
    hess[95] = t511;
    hess[96] = t229;
    hess[97] = t313;
    hess[98] = t376;
    hess[99] = t404;
    hess[100] = t420;
    hess[101] = t430;
    hess[102] = t466;
    hess[103] = t503;
    hess[104] = t119 * (-t442 - t512 - t513 + t514 + t516 - t521);
    hess[105] = t523;
    hess[106] = t525;
    hess[107] = t527;
    hess[108] = t232;
    hess[109] = t317;
    hess[110] = t380;
    hess[111] = t405;
    hess[112] = t421;
    hess[113] = t431;
    hess[114] = t472;
    hess[115] = t505;
    hess[116] = t523;
    hess[117] = t119 * (-t102 - t440 - t467);
    hess[118] = t528;
    hess[119] = t529;
    hess[120] = t235;
    hess[121] = t322;
    hess[122] = t385;
    hess[123] = t406;
    hess[124] = t422;
    hess[125] = t432;
    hess[126] = t476;
    hess[127] = t508;
    hess[128] = t525;
    hess[129] = t528;
    hess[130] = t119 * (-t104 - t490 + t506);
    hess[131] = t530;
    hess[132] = t238;
    hess[133] = t326;
    hess[134] = t388;
    hess[135] = t407;
    hess[136] = t423;
    hess[137] = t433;
    hess[138] = t479;
    hess[139] = t511;
    hess[140] = t527;
    hess[141] = t529;
    hess[142] = t530;
    hess[143] = t119 * (-t106 - t521 + t526);
}

// hess is (144×1) flattened in column-major order
void edge_edge_closest_point_hessian_b(
    double ea0_x,
    double ea0_y,
    double ea0_z,
    double ea1_x,
    double ea1_y,
    double ea1_z,
    double eb0_x,
    double eb0_y,
    double eb0_z,
    double eb1_x,
    double eb1_y,
    double eb1_z,
    double hess[144])
{
    const auto t0 = ea0_z * eb0_z;
    const auto t1 = ea1_z * eb1_z;
    const auto t2 = ea0_z * eb1_z;
    const auto t3 = ea1_z * eb0_z;
    const auto t4 = t0 + t1 - t2 - t3;
    const auto t5 = ea0_y * eb0_y;
    const auto t6 = ea1_y * eb1_y;
    const auto t7 = ea0_y * eb1_y;
    const auto t8 = ea1_y * eb0_y;
    const auto t9 = t5 + t6 - t7 - t8;
    const auto t10 = t4 + t9;
    const auto t11 = ea0_x * eb0_x;
    const auto t12 = ea1_x * eb1_x;
    const auto t13 = ea0_x * eb1_x;
    const auto t14 = ea1_x * eb0_x;
    const auto t15 = t11 + t12 - t13 - t14;
    const auto t16 = t10 + t15;
    const auto t17 = ea0_x - ea1_x;
    const auto t18 = ea0_x - eb0_x;
    const auto t19 = t17 * t18;
    const auto t20 = ea0_y - ea1_y;
    const auto t21 = ea0_y - eb0_y;
    const auto t22 = t20 * t21;
    const auto t23 = ea0_z - ea1_z;
    const auto t24 = ea0_z - eb0_z;
    const auto t25 = t23 * t24;
    const auto t26 = t22 + t25;
    const auto t27 = t19 + t26;
    const auto t28 = t16 * t27;
    const auto t29 = 2 * t11;
    const auto t30 = 2 * t13;
    const auto t31 = ea0_x * eb0_y;
    const auto t32 = ea1_x * eb1_y;
    const auto t33 = ea0_x * eb0_z;
    const auto t34 = ea1_x * eb1_z;
    const auto t35 = 2 * t5;
    const auto t36 = 2 * t7;
    const auto t37 = ea0_y * eb1_x;
    const auto t38 = ea1_y * eb0_x;
    const auto t39 = ea0_y * eb0_z;
    const auto t40 = ea1_y * eb1_z;
    const auto t41 = 2 * t0;
    const auto t42 = 2 * t2;
    const auto t43 = ea0_z * eb1_x;
    const auto t44 = ea1_z * eb0_x;
    const auto t45 = ea0_z * eb1_y;
    const auto t46 = ea1_z * eb0_y;
    const auto t47 = 2 * t14;
    const auto t48 = 2 * t12;
    const auto t49 = 2 * t8;
    const auto t50 = 2 * t6;
    const auto t51 = std::pow(ea0_x, 2);
    const auto t52 = std::pow(eb0_y, 2);
    const auto t53 = std::pow(eb0_z, 2);
    const auto t54 = std::pow(eb1_y, 2);
    const auto t55 = std::pow(eb1_z, 2);
    const auto t56 = std::pow(ea0_y, 2);
    const auto t57 = std::pow(eb0_x, 2);
    const auto t58 = std::pow(eb1_x, 2);
    const auto t59 = std::pow(ea0_z, 2);
    const auto t60 = std::pow(ea1_x, 2);
    const auto t61 = std::pow(ea1_y, 2);
    const auto t62 = std::pow(ea1_z, 2);
    const auto t63 = 2 * ea0_x;
    const auto t64 = ea1_x * t52;
    const auto t65 = ea1_x * t53;
    const auto t66 = ea1_x * t54;
    const auto t67 = ea1_x * t55;
    const auto t68 = 2 * eb0_y;
    const auto t69 = eb1_y * t51;
    const auto t70 = 2 * eb0_z;
    const auto t71 = eb1_z * t51;
    const auto t72 = 2 * ea0_y;
    const auto t73 = ea1_y * t57;
    const auto t74 = ea1_y * t53;
    const auto t75 = ea1_y * t58;
    const auto t76 = ea1_y * t55;
    const auto t77 = 2 * eb0_x;
    const auto t78 = eb1_x * t56;
    const auto t79 = eb1_z * t56;
    const auto t80 = 2 * ea0_z;
    const auto t81 = ea1_z * t57;
    const auto t82 = ea1_z * t52;
    const auto t83 = ea1_z * t58;
    const auto t84 = ea1_z * t54;
    const auto t85 = eb1_x * t59;
    const auto t86 = eb1_y * t59;
    const auto t87 = eb1_y * t60;
    const auto t88 = eb1_z * t60;
    const auto t89 = eb1_x * t61;
    const auto t90 = eb1_z * t61;
    const auto t91 = eb1_x * t62;
    const auto t92 = eb1_y * t62;
    const auto t93 = -t0 * t29 + t0 * t30 - t0 * t35 + t0 * t36 - t1 * t29
        + t1 * t30 - t1 * t35 + t1 * t36 + t1 * t47 - t1 * t48 + t1 * t49
        - t1 * t50 - t12 * t35 + t12 * t36 - t12 * t41 + t12 * t42 + t14 * t35
        - t14 * t36 + t14 * t41 - t14 * t42 + t2 * t29 - t2 * t30 + t2 * t35
        - t2 * t36 + t29 * t3 - t29 * t5 - t29 * t6 + t29 * t7 + t29 * t8
        - t3 * t30 + t3 * t35 - t3 * t36 - t3 * t47 + t3 * t48 - t3 * t49
        + t3 * t50 + t30 * t5 + t30 * t6 - t30 * t7 - t30 * t8 + 4 * t31 * t32
        + 4 * t33 * t34 + 4 * t37 * t38 + 4 * t39 * t40 - t41 * t6 + t41 * t8
        + t42 * t6 - t42 * t8 + 4 * t43 * t44 + 4 * t45 * t46 + t47 * t6
        - t47 * t8 - t48 * t6 + t48 * t8 + t51 * t52 + t51 * t53 + t51 * t54
        + t51 * t55 + t52 * t59 + t52 * t60 + t52 * t62 + t53 * t56 + t53 * t60
        + t53 * t61 + t54 * t59 + t54 * t60 + t54 * t62 + t55 * t56 + t55 * t60
        + t55 * t61 + t56 * t57 + t56 * t58 + t57 * t59 + t57 * t61 + t57 * t62
        + t58 * t59 + t58 * t61 + t58 * t62 - t63 * t64 - t63 * t65 - t63 * t66
        - t63 * t67 - t68 * t69 - t68 * t86 - t68 * t87 - t68 * t92 - t70 * t71
        - t70 * t79 - t70 * t88 - t70 * t90 - t72 * t73 - t72 * t74 - t72 * t75
        - t72 * t76 - t77 * t78 - t77 * t85 - t77 * t89 - t77 * t91 - t80 * t81
        - t80 * t82 - t80 * t83 - t80 * t84;
    const auto t94 = 1.0 / t93;
    const auto t95 = -eb1_z * t70 + t53 + t55;
    const auto t96 = -eb1_y * t68 + t52 + t54;
    const auto t97 = t95 + t96;
    const auto t98 = ea1_z * t80;
    const auto t99 = t59 + t62 - t98;
    const auto t100 = ea1_y * t72;
    const auto t101 = -t100 + t56 + t61;
    const auto t102 = t101 + t99;
    const auto t103 = ea1_x * t63;
    const auto t104 = -t103 + t51 + t60;
    const auto t105 = t102 + t104;
    const auto t106 = eb0_x - eb1_x;
    const auto t107 = t106 * t18;
    const auto t108 = eb0_y - eb1_y;
    const auto t109 = t108 * t21;
    const auto t110 = eb0_z - eb1_z;
    const auto t111 = t110 * t24;
    const auto t112 = t109 + t111;
    const auto t113 = t107 + t112;
    const auto t114 = eb1_y * t63;
    const auto t115 = eb1_z * t63;
    const auto t116 = ea0_x * t52 + ea0_x * t53 + ea0_x * t54 + ea0_x * t55
        - eb0_x * t0 - eb0_x * t1 + eb0_x * t2 + eb0_x * t3 - eb0_x * t5
        - eb0_x * t6 + eb0_x * t7 + eb0_x * t8 - eb0_y * t114 - eb0_z * t115
        + eb1_x * t0 + eb1_x * t1 - eb1_x * t2 - eb1_x * t3 + eb1_x * t5
        + eb1_x * t6 - eb1_x * t7 - eb1_x * t8 + t32 * t68 + t34 * t70 - t64
        - t65 - t66 - t67;
    const auto t117 = std::pow(t116, 2);
    const auto t118 = std::pow(t93, -2);
    const auto t119 = 2 * t94;
    const auto t120 = t116 * t119;
    const auto t121 = t106 * t120;
    const auto t122 = t117 * t118;
    const auto t123 = t105 * t113;
    const auto t124 = 4 * t123;
    const auto t125 = -t105 * t113 * t94 * t97 - 4 * t113 * t116 * t17 * t94
        - 4 * t117 * t118 * t16 * t27 + t121 * t27 + t122 * t124
        + t28 * t94 * t97;
    const auto t126 =
        -t0 - t1 - t11 - t12 + t13 + t14 + t2 + t3 - t5 - t6 + t7 + t8;
    const auto t127 = t113 + t126;
    const auto t128 = ea1_x + eb0_x - t63;
    const auto t129 = t106 * t17;
    const auto t130 = t128 * t16;
    const auto t131 = -t105 * t121 + t106 * t128 - t120 * t130 + 2 * t129;
    const auto t132 = ea1_y + eb0_y - t72;
    const auto t133 = t106 * t132;
    const auto t134 = t108 * t17;
    const auto t135 = 2 * t134;
    const auto t136 = -t17 * t18 - t20 * t21 - t23 * t24;
    const auto t137 = t108 * t120;
    const auto t138 = t136 * t137;
    const auto t139 = t105 * t137;
    const auto t140 = t120 * t16;
    const auto t141 = t132 * t140;
    const auto t142 = eb0_x * t72;
    const auto t143 = eb1_z * t72;
    const auto t144 = ea1_y * eb1_x;
    const auto t145 = -ea0_y * t53 - ea0_y * t55 - ea0_y * t57 - ea0_y * t58
        + eb0_y * t0 + eb0_y * t1 + eb0_y * t11 + eb0_y * t12 - eb0_y * t13
        - eb0_y * t14 - eb0_y * t2 - eb0_y * t3 + eb0_z * t143 + eb1_x * t142
        - eb1_y * t0 - eb1_y * t1 - eb1_y * t11 - eb1_y * t12 + eb1_y * t13
        + eb1_y * t14 + eb1_y * t2 + eb1_y * t3 - t144 * t77 - t40 * t70 + t73
        + t74 + t75 + t76;
    const auto t146 = t119 * t145;
    const auto t147 = t106 * t146;
    const auto t148 = t136 * t147;
    const auto t149 = t136 * t16;
    const auto t150 =
        t119 * (eb0_x * eb0_y - eb0_x * eb1_y - eb0_y * eb1_x + eb1_x * eb1_y);
    const auto t151 = t149 * t150;
    const auto t152 = 8 * t116;
    const auto t153 = t118 * t152;
    const auto t154 = t145 * t153;
    const auto t155 = t149 * t154;
    const auto t156 = t106 * t20;
    const auto t157 = 4 * t116;
    const auto t158 = t113 * t94;
    const auto t159 = t157 * t158;
    const auto t160 = t159 * t20;
    const auto t161 = t123 * t150;
    const auto t162 = t158 * t17;
    const auto t163 = 4 * t162;
    const auto t164 = t145 * t163;
    const auto t165 = t123 * t154;
    const auto t166 = t105 * t147 + t108 * t128 + t130 * t146 + 2 * t156 - t160
        + t161 + t164 - t165;
    const auto t167 =
        t94 * (-t133 - t135 + t138 + t139 + t141 - t148 - t151 + t155 - t166);
    const auto t168 = ea1_z + eb0_z - t80;
    const auto t169 = t106 * t168;
    const auto t170 = t110 * t17;
    const auto t171 = 2 * t170;
    const auto t172 = t110 * t120;
    const auto t173 = t136 * t172;
    const auto t174 = t105 * t172;
    const auto t175 = t140 * t168;
    const auto t176 = eb0_x * t80;
    const auto t177 = eb0_y * t80;
    const auto t178 = ea1_z * eb1_x;
    const auto t179 = ea1_z * eb1_y;
    const auto t180 = -ea0_z * t52 - ea0_z * t54 - ea0_z * t57 - ea0_z * t58
        + eb0_z * t11 + eb0_z * t12 - eb0_z * t13 - eb0_z * t14 + eb0_z * t5
        + eb0_z * t6 - eb0_z * t7 - eb0_z * t8 + eb1_x * t176 + eb1_y * t177
        - eb1_z * t11 - eb1_z * t12 + eb1_z * t13 + eb1_z * t14 - eb1_z * t5
        - eb1_z * t6 + eb1_z * t7 + eb1_z * t8 - t178 * t77 - t179 * t68 + t81
        + t82 + t83 + t84;
    const auto t181 = t119 * t180;
    const auto t182 = t106 * t181;
    const auto t183 = t136 * t182;
    const auto t184 =
        t119 * (eb0_x * eb0_z - eb0_x * eb1_z - eb0_z * eb1_x + eb1_x * eb1_z);
    const auto t185 = t149 * t184;
    const auto t186 = t153 * t180;
    const auto t187 = t149 * t186;
    const auto t188 = t106 * t23;
    const auto t189 = t159 * t23;
    const auto t190 = t123 * t184;
    const auto t191 = t163 * t180;
    const auto t192 = t123 * t186;
    const auto t193 = t105 * t182 + t110 * t128 + t130 * t181 + 2 * t188 - t189
        + t190 + t191 - t192;
    const auto t194 =
        t94 * (-t169 - t171 + t173 + t174 + t175 - t183 - t185 + t187 - t193);
    const auto t195 = t140 * t18;
    const auto t196 = 2 * t109;
    const auto t197 = t136 * t94;
    const auto t198 = t119 * t97;
    const auto t199 = 8 * t122;
    const auto t200 = 2 * t111 + t126;
    const auto t201 = t94
        * (-t106 * t157 * t197 + t107 - t123 * t198 + t123 * t199 + t131
           - t149 * t198 + t149 * t199 - t152 * t162 + t195 + t196 + t200);
    const auto t202 = t106 * t21;
    const auto t203 = t147 * t27;
    const auto t204 = t150 * t28;
    const auto t205 = t140 * t21;
    const auto t206 = t137 * t27;
    const auto t207 = t154 * t28;
    const auto t208 = t94 * (t166 - t202 - t203 - t204 + t205 + t206 + t207);
    const auto t209 = t106 * t24;
    const auto t210 = t182 * t27;
    const auto t211 = t184 * t28;
    const auto t212 = t140 * t24;
    const auto t213 = t172 * t27;
    const auto t214 = t186 * t28;
    const auto t215 = t94 * (t193 - t209 - t210 - t211 + t212 + t213 + t214);
    const auto t216 = ea0_x * t0 + ea0_x * t1 - ea0_x * t2 - ea0_x * t3
        + ea0_x * t5 + ea0_x * t6 - ea0_x * t7 - ea0_x * t8 - ea1_x * t0
        - ea1_x * t1 + ea1_x * t2 + ea1_x * t3 - ea1_x * t5 - ea1_x * t6
        + ea1_x * t7 + ea1_x * t8 - eb0_x * t56 - eb0_x * t59 - eb0_x * t61
        - eb0_x * t62 - t144 * t72 - t178 * t80 + t38 * t72 + t44 * t80 + t78
        + t85 + t89 + t91;
    const auto t217 = t105 * t216;
    const auto t218 = t106 * t119;
    const auto t219 = t119 * t130;
    const auto t220 = t10 * t119;
    const auto t221 = t123 * t220;
    const auto t222 = t163 * t216;
    const auto t223 = t153 * t216;
    const auto t224 = t123 * t223;
    const auto t225 = t221 + t222 - t224;
    const auto t226 = t128 * t17 + t216 * t219 + t217 * t218 + t225;
    const auto t227 = ea0_x + eb1_x - t77;
    const auto t228 = 2 * t17;
    const auto t229 = t120 * t136;
    const auto t230 = t105 * t120;
    const auto t231 = t136 * t216;
    const auto t232 = t129 + t136 - t140 * t17 + t149 * t220 - t149 * t223 + t16
        - t17 * t229 + t218 * t231 + t227 * t228 - t227 * t230;
    const auto t233 = t94 * (-t105 - t226 - t232);
    const auto t234 = t128 * t20;
    const auto t235 = ea0_y + eb1_y - t68;
    const auto t236 = t228 * t235;
    const auto t237 = t140 * t20;
    const auto t238 = ea1_x * eb0_y;
    const auto t239 = -ea0_y * t0 - ea0_y * t1 - ea0_y * t11 - ea0_y * t12
        + ea0_y * t13 + ea0_y * t14 + ea0_y * t2 + ea0_y * t3 + ea1_y * t0
        + ea1_y * t1 + ea1_y * t11 + ea1_y * t12 - ea1_y * t13 - ea1_y * t14
        - ea1_y * t2 - ea1_y * t3 + eb0_y * t51 + eb0_y * t59 + eb0_y * t60
        + eb0_y * t62 + t179 * t80 - t238 * t63 + t32 * t63 - t46 * t80 - t69
        - t86 - t87 - t92;
    const auto t240 = t105 * t239;
    const auto t241 = t218 * t240;
    const auto t242 = t230 * t235;
    const auto t243 = t219 * t239;
    const auto t244 = t17 * t239;
    const auto t245 = 4 * t158;
    const auto t246 = t244 * t245;
    const auto t247 = t119
        * (-ea0_y * eb0_x - ea1_x * t68 + eb0_y * t63 - t114 - t144 + 2 * t32
           + t37 + t38);
    const auto t248 = t123 * t247;
    const auto t249 = t120 * t27;
    const auto t250 = t239 * t27;
    const auto t251 = t123 * t153 * t239;
    const auto t252 = -8 * t116 * t118 * t16 * t239 * t27 + t20 * t249
        + t218 * t250 - t246 + t247 * t28 - t248 + t251;
    const auto t253 =
        t94 * (-t156 - t234 - t236 + t237 + t241 + t242 + t243 - t252);
    const auto t254 = t128 * t23;
    const auto t255 = ea0_z + eb1_z - t70;
    const auto t256 = t228 * t255;
    const auto t257 = t140 * t23;
    const auto t258 = ea1_x * eb0_z;
    const auto t259 = ea1_y * eb0_z;
    const auto t260 = -ea0_z * t11 - ea0_z * t12 + ea0_z * t13 + ea0_z * t14
        - ea0_z * t5 - ea0_z * t6 + ea0_z * t7 + ea0_z * t8 + ea1_z * t11
        + ea1_z * t12 - ea1_z * t13 - ea1_z * t14 + ea1_z * t5 + ea1_z * t6
        - ea1_z * t7 - ea1_z * t8 + eb0_z * t51 + eb0_z * t56 + eb0_z * t60
        + eb0_z * t61 - t258 * t63 - t259 * t72 + t34 * t63 + t40 * t72 - t71
        - t79 - t88 - t90;
    const auto t261 = t105 * t260;
    const auto t262 = t218 * t261;
    const auto t263 = t230 * t255;
    const auto t264 = t219 * t260;
    const auto t265 = t17 * t260;
    const auto t266 = t245 * t265;
    const auto t267 = t119
        * (-ea0_z * eb0_x - ea1_x * t70 + eb0_z * t63 - t115 - t178 + 2 * t34
           + t43 + t44);
    const auto t268 = t123 * t267;
    const auto t269 = t260 * t27;
    const auto t270 = t123 * t153 * t260;
    const auto t271 = -8 * t116 * t118 * t16 * t260 * t27 + t218 * t269
        + t23 * t249 - t266 + t267 * t28 - t268 + t270;
    const auto t272 =
        t94 * (-t188 - t254 - t256 + t257 + t262 + t263 + t264 - t271);
    const auto t273 = -t22;
    const auto t274 = t105 * t18;
    const auto t275 = t120 * t274;
    const auto t276 = t216 * t27;
    const auto t277 = t218 * t276;
    const auto t278 = t220 * t28;
    const auto t279 = t17 * t27;
    const auto t280 = t120 * t279;
    const auto t281 = t16 * t276;
    const auto t282 = t153 * t281;
    const auto t283 = -t25;
    const auto t284 = t105 + t283;
    const auto t285 =
        t94 * (t19 + t226 + t273 - t275 - t277 - t278 + t280 + t282 + t284);
    const auto t286 = t17 * t21;
    const auto t287 = 2 * t286;
    const auto t288 = t21 * t230;
    const auto t289 = t136 * t218;
    const auto t290 = -8 * t116 * t118 * t136 * t16 * t239 + t149 * t247
        + t20 * t229 + t239 * t289 + t246 + t248 - t251;
    const auto t291 = t94 * (t234 - t241 - t243 + t287 - t288 - t290);
    const auto t292 = t17 * t24;
    const auto t293 = 2 * t292;
    const auto t294 = t230 * t24;
    const auto t295 = -8 * t116 * t118 * t136 * t16 * t260 + t149 * t267
        + t229 * t23 + t260 * t289 + t266 + t268 - t270;
    const auto t296 = t94 * (t254 - t262 - t264 + t293 - t294 - t295);
    const auto t297 = -eb1_x * t77 + t57 + t58;
    const auto t298 = t297 + t95;
    const auto t299 = t108 * t146;
    const auto t300 = std::pow(t145, 2);
    const auto t301 = t118 * t300;
    const auto t302 = t149 * t301;
    const auto t303 = t20 * t245;
    const auto t304 = -t105 * t113 * t298 * t94 + t124 * t301 + t145 * t303;
    const auto t305 = t108 * t20;
    const auto t306 = t132 * t16;
    const auto t307 = t105 * t299 + t108 * t132 + t146 * t306 + 2 * t305;
    const auto t308 =
        t119 * (eb0_y * eb0_z - eb0_y * eb1_z - eb0_z * eb1_y + eb1_y * eb1_z);
    const auto t309 = t123 * t308;
    const auto t310 = t180 * t303;
    const auto t311 = t23 * t245;
    const auto t312 = t145 * t311;
    const auto t313 = 8 * t123;
    const auto t314 = t118 * t145;
    const auto t315 = t180 * t314;
    const auto t316 = t313 * t315;
    const auto t317 = t108 * t181;
    const auto t318 = t110 * t146;
    const auto t319 = 8 * t149;
    const auto t320 = t108 * t23;
    const auto t321 = t105 * t317 + t110 * t132 + t181 * t306 + 2 * t320;
    const auto t322 = t110 * t20;
    const auto t323 = t146 * t16;
    const auto t324 = t105 * t318 + t108 * t168 + t168 * t323 + 2 * t322;
    const auto t325 = -t94
        * (t136 * t317 + t136 * t318 + t149 * t308 + t309 + t310 + t312
           + t315 * t319 + t316 + t321 + t324);
    const auto t326 = t108 * t18 + t160 - t161 - t164 + t165 + t18 * t323;
    const auto t327 =
        t94 * (t133 + t135 - t138 - t139 - t141 + t148 + t151 - t155 - t326);
    const auto t328 = -t21 * t323;
    const auto t329 = 2 * t107;
    const auto t330 = t119 * t298;
    const auto t331 = 4 * t197;
    const auto t332 = 8 * t158;
    const auto t333 = t94
        * (t108 * t145 * t331 + t109 - t123 * t330 + t145 * t20 * t332
           - t149 * t330 + t200 + t301 * t313 + 8 * t302 + t307 + t328 + t329);
    const auto t334 = t108 * t24;
    const auto t335 = t24 * t323;
    const auto t336 = t27 * t317;
    const auto t337 = t27 * t318;
    const auto t338 = t28 * t308;
    const auto t339 = 8 * t28;
    const auto t340 = t315 * t339;
    const auto t341 = t309 + t310 + t312 + t316 - t336 - t337 - t338 - t340;
    const auto t342 = t94 * (t321 - t334 - t335 + t341);
    const auto t343 = t136 * t146;
    const auto t344 = t108 * t119;
    const auto t345 = t119
        * (-ea0_x * eb1_y + ea1_y * t77 + eb1_x * t72 - t142 - 2 * t144 - t238
           + t31 + t32);
    const auto t346 = t216 * t314;
    const auto t347 = t149 * t345 + t17 * t343 + t231 * t344 + t319 * t346;
    const auto t348 = t119 * t306;
    const auto t349 = t132 * t17 + t216 * t348 + t217 * t344;
    const auto t350 = t123 * t345;
    const auto t351 = t20 * t216;
    const auto t352 = t245 * t351;
    const auto t353 = t313 * t346;
    const auto t354 = 2 * t20;
    const auto t355 = t105 * t146;
    const auto t356 =
        t134 + t17 * t323 + t227 * t354 + t227 * t355 + t350 + t352 + t353;
    const auto t357 = -t94 * (t347 + t349 + t356);
    const auto t358 = t119 * (t15 + t4);
    const auto t359 = -2 * t108 * t136 * t239 * t94
        - 8 * t118 * t136 * t145 * t16 * t239 + t136 + t149 * t358 + t20 * t343;
    const auto t360 = -t239 * t303;
    const auto t361 = t123 * t358;
    const auto t362 = t239 * t314;
    const auto t363 = -t313 * t362;
    const auto t364 =
        t132 * t20 - t239 * t348 - t240 * t344 + t360 + t361 + t363;
    const auto t365 = t16 + t20 * t323 + t235 * t354 + t235 * t355 + t305;
    const auto t366 = t94 * (-t105 - t359 - t364 - t365);
    const auto t367 = t132 * t23;
    const auto t368 = t255 * t354;
    const auto t369 = t23 * t323;
    const auto t370 = t255 * t355;
    const auto t371 = t261 * t344;
    const auto t372 = t260 * t348;
    const auto t373 = t119
        * (-ea0_z * eb0_y - ea1_y * t70 + eb0_z * t72 - t143 - t179 + 2 * t40
           + t45 + t46);
    const auto t374 = t123 * t373;
    const auto t375 = t20 * t260;
    const auto t376 = t245 * t375;
    const auto t377 = t260 * t314;
    const auto t378 = t313 * t377;
    const auto t379 = t136 * t260 * t344 + t149 * t373 - t23 * t343
        + t319 * t377 + t374 + t376 + t378;
    const auto t380 =
        t94 * (-t320 - t367 - t368 - t369 - t370 + t371 + t372 + t379);
    const auto t381 = t18 * t20;
    const auto t382 = t146 * t274 + t350 + t352 + t353 + 2 * t381;
    const auto t383 = 8 * t281;
    const auto t384 = -t146 * t279 - t276 * t344 - t28 * t345 - t314 * t383;
    const auto t385 = t94 * (t349 + t382 + t384);
    const auto t386 = t21 * t355 + t22;
    const auto t387 = -t19;
    const auto t388 = t146 * t27;
    const auto t389 =
        -t20 * t388 + t250 * t344 - t28 * t358 + t339 * t362 + t387;
    const auto t390 = t94 * (t284 + t364 + t386 + t389);
    const auto t391 = t20 * t24;
    const auto t392 = 2 * t391;
    const auto t393 = t24 * t355;
    const auto t394 = -t23 * t388 + t269 * t344 + t28 * t373 + t339 * t377
        - t374 - t376 - t378;
    const auto t395 = t94 * (t367 - t371 - t372 + t392 + t393 + t394);
    const auto t396 = t297 + t96;
    const auto t397 = t110 * t181;
    const auto t398 = std::pow(t180, 2);
    const auto t399 = t118 * t398;
    const auto t400 = t149 * t399;
    const auto t401 = -t105 * t113 * t396 * t94 + t124 * t399 + t180 * t311;
    const auto t402 = t110 * t23;
    const auto t403 = t16 * t181;
    const auto t404 = t105 * t397 + t110 * t168 + t168 * t403 + 2 * t402;
    const auto t405 = t110 * t18 + t18 * t403 + t189 - t190 - t191 + t192;
    const auto t406 =
        t94 * (t169 + t171 - t173 - t174 - t175 + t183 + t185 - t187 - t405);
    const auto t407 = t110 * t21;
    const auto t408 = t21 * t403;
    const auto t409 = t94 * (t324 + t341 - t407 - t408);
    const auto t410 = -t24 * t403;
    const auto t411 = t119 * t396;
    const auto t412 = t94
        * (t110 * t180 * t331 + t111 - t123 * t411 + t126 - t149 * t411
           + t180 * t23 * t332 + t196 + t313 * t399 + t329 + 8 * t400 + t404
           + t410);
    const auto t413 = t136 * t181;
    const auto t414 = t110 * t119;
    const auto t415 = t119
        * (-ea0_x * eb1_z + ea1_z * t77 + eb1_x * t80 - t176 - 2 * t178 - t258
           + t33 + t34);
    const auto t416 = t118 * t180;
    const auto t417 = t216 * t416;
    const auto t418 = t149 * t415 + t17 * t413 + t231 * t414 + t319 * t417;
    const auto t419 = t119 * t16;
    const auto t420 = t168 * t419;
    const auto t421 = t168 * t17 + t216 * t420 + t217 * t414;
    const auto t422 = t123 * t415;
    const auto t423 = t216 * t23;
    const auto t424 = t245 * t423;
    const auto t425 = t313 * t417;
    const auto t426 = 2 * t23;
    const auto t427 = t105 * t181;
    const auto t428 =
        t17 * t403 + t170 + t227 * t426 + t227 * t427 + t422 + t424 + t425;
    const auto t429 = -t94 * (t418 + t421 + t428);
    const auto t430 = t119
        * (-ea0_y * eb1_z + ea1_z * t68 + eb1_y * t80 - t177 - 2 * t179 - t259
           + t39 + t40);
    const auto t431 = -2 * t110 * t136 * t239 * t94
        - 8 * t118 * t136 * t16 * t180 * t239 + t149 * t430 + t20 * t413;
    const auto t432 = t168 * t20 - t239 * t420 - t240 * t414;
    const auto t433 = t23 * t239;
    const auto t434 = -t245 * t433;
    const auto t435 = t123 * t430;
    const auto t436 = t239 * t416;
    const auto t437 = -t313 * t436;
    const auto t438 =
        t20 * t403 + t235 * t426 + t235 * t427 + t322 + t434 + t435 + t437;
    const auto t439 = t94 * (-t431 - t432 - t438);
    const auto t440 = t119 * (t15 + t9);
    const auto t441 = -2 * t110 * t136 * t260 * t94
        - 8 * t118 * t136 * t16 * t180 * t260 + t136 + t149 * t440 + t23 * t413;
    const auto t442 = -t260 * t311;
    const auto t443 = t123 * t440;
    const auto t444 = t260 * t416;
    const auto t445 = -t313 * t444;
    const auto t446 =
        t105 + t168 * t23 - t260 * t420 - t261 * t414 + t442 + t443 + t445;
    const auto t447 = t16 + t23 * t403 + t255 * t426 + t255 * t427 + t402;
    const auto t448 = t94 * (-t441 - t446 - t447);
    const auto t449 = t18 * t23;
    const auto t450 = t181 * t274 + t422 + t424 + t425 + 2 * t449;
    const auto t451 = -t181 * t279 - t276 * t414 - t28 * t415 - t383 * t416;
    const auto t452 = t94 * (t421 + t450 + t451);
    const auto t453 = t21 * t23;
    const auto t454 = t21 * t427 + t434 + t435 + t437 + 2 * t453;
    const auto t455 = t181 * t27;
    const auto t456 = -t20 * t455 + t250 * t414 - t28 * t430 + t339 * t436;
    const auto t457 = t94 * (t432 + t454 + t456);
    const auto t458 = t24 * t427 + t25;
    const auto t459 =
        -t23 * t455 + t269 * t414 + t273 - t28 * t440 + t339 * t444 + t387;
    const auto t460 = t94 * (t446 + t458 + t459);
    const auto t461 = t94 * (t202 + t203 + t204 - t205 - t206 - t207 + t326);
    const auto t462 = t94 * (t209 + t210 + t211 - t212 - t213 - t214 + t405);
    const auto t463 = t18 * t419;
    const auto t464 = t216 * t463;
    const auto t465 = t94 * (t225 + t232 + t387 - t464);
    const auto t466 = t239 * t463;
    const auto t467 = t94 * (t156 + t236 - t237 - t242 - t290 - t381 + t466);
    const auto t468 = t260 * t463;
    const auto t469 = t94 * (t188 + t256 - t257 - t263 - t295 - t449 + t468);
    const auto t470 = t94
        * (-t221 - t222 + t224 + t26 + t275 + t277 + t278 - t280 - t282 + t464);
    const auto t471 = t94 * (-t252 - t287 + t288 + t381 - t466);
    const auto t472 = t94 * (-t271 - t293 + t294 + t449 - t468);
    const auto t473 = t94
        * (-t309 - t310 - t312 - t316 + t334 + t335 + t336 + t337 + t338 + t340
           + t407 + t408);
    const auto t474 = t21 * t419;
    const auto t475 = -t216 * t474 - t286;
    const auto t476 = t94 * (t356 + t384 + t475);
    const auto t477 = t239 * t474 + t360 + t361 + t363;
    const auto t478 = t94 * (-2 * t22 + t283 + t365 + t389 + t477);
    const auto t479 = t260 * t474;
    const auto t480 = t94 * (t320 + t368 + t369 + t370 + t394 - t453 + t479);
    const auto t481 = t94 * (-t347 - t382 - t475);
    const auto t482 = t94 * (-t359 - t386 - t477);
    const auto t483 = t94 * (t379 - t392 - t393 + t453 - t479);
    const auto t484 = t24 * t419;
    const auto t485 = -t216 * t484 - t292;
    const auto t486 = t94 * (t428 + t451 + t485);
    const auto t487 = t239 * t484 - t391;
    const auto t488 = t94 * (t438 + t456 + t487);
    const auto t489 = t260 * t484 + t442 + t443 + t445;
    const auto t490 = t94 * (-2 * t25 + t447 + t459 + t489);
    const auto t491 = t94 * (-t418 - t450 - t485);
    const auto t492 = t94 * (-t431 - t454 - t487);
    const auto t493 = t94 * (-t441 - t458 - t489);
    const auto t494 = std::pow(t17, 2);
    const auto t495 = t119 * t227;
    const auto t496 = t217 * t495;
    const auto t497 = t17 * t216 * t419;
    const auto t498 = std::pow(t216, 2);
    const auto t499 = t118 * t498;
    const auto t500 = t17 * t20;
    const auto t501 =
        ea0_x * ea0_y - ea0_x * ea1_y - ea0_y * ea1_x + ea1_x * ea1_y;
    const auto t502 = t123 * t501;
    const auto t503 = t149 * t94;
    const auto t504 = t217 * t94;
    const auto t505 = t16 * t94;
    const auto t506 = t119
        * (4 * t105 * t113 * t118 * t216 * t239 + t105 * t227 * t239 * t94
           + 4 * t118 * t136 * t16 * t216 * t239 + t136 * t17 * t239 * t94
           + t16 * t17 * t239 * t94 - t197 * t351 - t235 * t504 - t351 * t505
           - t500 - t501 * t503 - t502 * t94);
    const auto t507 = t17 * t23;
    const auto t508 =
        ea0_x * ea0_z - ea0_x * ea1_z - ea0_z * ea1_x + ea1_x * ea1_z;
    const auto t509 = t123 * t508;
    const auto t510 = t119
        * (4 * t105 * t113 * t118 * t216 * t260 + t105 * t227 * t260 * t94
           + 4 * t118 * t136 * t16 * t216 * t260 + t136 * t17 * t260 * t94
           + t16 * t17 * t260 * t94 - t197 * t423 - t255 * t504 - t423 * t505
           - t503 * t508 - t507 - t509 * t94);
    const auto t511 = 4 * t276;
    const auto t512 = t94
        * (-t102 * t119 * t123 + 2 * t102 * t16 * t27 * t94
           + 8 * t105 * t113 * t118 * t498 + 2 * t105 * t18 * t216 * t94 - t105
           - t17 * t511 * t94 - t339 * t499 + t494 + t496 + t497);
    const auto t513 = t21 * t217;
    const auto t514 = t239 * t279;
    const auto t515 = t118 * t216;
    const auto t516 = t118 * t383;
    const auto t517 = -t119 * t20 * t276 - t119 * t28 * t501 + t119 * t502
        + t119 * t514 - t239 * t313 * t515 + t239 * t516 + t500;
    const auto t518 = t94 * (t119 * t513 - t240 * t495 - t244 * t419 + t517);
    const auto t519 = t217 * t24;
    const auto t520 = t260 * t279;
    const auto t521 = -t119 * t23 * t276 - t119 * t28 * t508 + t119 * t509
        + t119 * t520 - t260 * t313 * t515 + t260 * t516 + t507;
    const auto t522 = t94 * (t119 * t519 - t261 * t495 - t265 * t419 + t521);
    const auto t523 = t104 + t99;
    const auto t524 = t123 * t523;
    const auto t525 = t28 * t523;
    const auto t526 = std::pow(t239, 2);
    const auto t527 = t118 * t526;
    const auto t528 = 4 * t28;
    const auto t529 = t20 * t239;
    const auto t530 = t119 * t235;
    const auto t531 = t105 - std::pow(t20, 2) + t240 * t530 + t419 * t529;
    const auto t532 = t20 * t23;
    const auto t533 =
        ea0_y * ea0_z - ea0_y * ea1_z - ea0_z * ea1_y + ea1_y * ea1_z;
    const auto t534 = t123 * t533;
    const auto t535 = t235 * t261;
    const auto t536 = t240 * t255;
    const auto t537 = t27 * t375;
    const auto t538 = t27 * t433;
    const auto t539 = t239 * t260;
    const auto t540 = t118 * t539;
    const auto t541 = t119
        * (-t124 * t540 + t28 * t533 * t94 + t375 * t505 + t433 * t505
           + t528 * t540 - t532 - t534 * t94 + t535 * t94 + t536 * t94
           - t537 * t94 - t538 * t94);
    const auto t542 =
        t94 * (-t119 * t239 * t274 + t217 * t530 + t351 * t419 + t517);
    const auto t543 = t119 * t149;
    const auto t544 = t94
        * (8 * t105 * t113 * t118 * t526 + 8 * t118 * t136 * t16 * t526
           - t119 * t21 * t240 - t119 * t524 - t331 * t529 - t523 * t543
           - t531);
    const auto t545 = t119 * t136;
    const auto t546 = t119 * t534 + t313 * t540 + t319 * t540 - t375 * t545
        - t433 * t545 + t532 + t533 * t543;
    const auto t547 =
        t94 * (-t119 * t24 * t240 - t119 * t535 - t375 * t419 + t546);
    const auto t548 = t101 + t104;
    const auto t549 = t123 * t548;
    const auto t550 = t28 * t548;
    const auto t551 = std::pow(t260, 2);
    const auto t552 = t118 * t551;
    const auto t553 = t23 * t260;
    const auto t554 = t119 * t255;
    const auto t555 = t105 - std::pow(t23, 2) + t261 * t554 + t419 * t553;
    const auto t556 =
        t94 * (-t119 * t260 * t274 + t217 * t554 + t419 * t423 + t521);
    const auto t557 =
        t94 * (-t119 * t21 * t261 - t119 * t536 - t419 * t433 + t546);
    const auto t558 = t94
        * (8 * t105 * t113 * t118 * t551 + 8 * t118 * t136 * t16 * t551
           - t119 * t24 * t261 - t119 * t549 - t331 * t553 - t543 * t548
           - t555);
    const auto t559 = t124 * t94;
    const auto t560 = 2 * t118;
    const auto t561 = t505 * t511;
    const auto t562 = t560
        * (4 * t105 * t113 * t216 * t239 * t94 + t105 * t18 * t239
           + t16 * t27 * t501 + t20 * t216 * t27 - t239 * t561 - t502 - t513
           - t514);
    const auto t563 = t560
        * (4 * t105 * t113 * t216 * t260 * t94 + t105 * t18 * t260
           + t16 * t27 * t508 + t216 * t23 * t27 - t260 * t561 - t509 - t519
           - t520);
    const auto t564 = t560
        * (t105 * t21 * t260 + t105 * t239 * t24
           + 4 * t16 * t239 * t260 * t27 * t94 + t16 * t27 * t533 - t534 - t537
           - t538 - t539 * t559);
    hess[0] = t119 * (-t125 - t127 - t131);
    hess[1] = t167;
    hess[2] = t194;
    hess[3] = t201;
    hess[4] = t208;
    hess[5] = t215;
    hess[6] = t233;
    hess[7] = t253;
    hess[8] = t272;
    hess[9] = t285;
    hess[10] = t291;
    hess[11] = t296;
    hess[12] = t167;
    hess[13] = t119
        * (-t127 + t136 * t16 * t298 * t94 - t136 * t299 - 4 * t302 - t304
           - t307);
    hess[14] = t325;
    hess[15] = t327;
    hess[16] = t333;
    hess[17] = t342;
    hess[18] = t357;
    hess[19] = t366;
    hess[20] = t380;
    hess[21] = t385;
    hess[22] = t390;
    hess[23] = t395;
    hess[24] = t194;
    hess[25] = t325;
    hess[26] = t119
        * (-t127 + t136 * t16 * t396 * t94 - t136 * t397 - 4 * t400 - t401
           - t404);
    hess[27] = t406;
    hess[28] = t409;
    hess[29] = t412;
    hess[30] = t429;
    hess[31] = t439;
    hess[32] = t448;
    hess[33] = t452;
    hess[34] = t457;
    hess[35] = t460;
    hess[36] = t201;
    hess[37] = t327;
    hess[38] = t406;
    hess[39] = t119 * (-t112 - t125 - t195);
    hess[40] = t461;
    hess[41] = t462;
    hess[42] = t465;
    hess[43] = t467;
    hess[44] = t469;
    hess[45] = t470;
    hess[46] = t471;
    hess[47] = t472;
    hess[48] = t208;
    hess[49] = t333;
    hess[50] = t409;
    hess[51] = t461;
    hess[52] = t119
        * (-t107 + 2 * t108 * t145 * t27 * t94 - t111
           + 4 * t118 * t16 * t27 * t300 - t28 * t298 * t94 - t304 - t328);
    hess[53] = t473;
    hess[54] = t476;
    hess[55] = t478;
    hess[56] = t480;
    hess[57] = t481;
    hess[58] = t482;
    hess[59] = t483;
    hess[60] = t215;
    hess[61] = t342;
    hess[62] = t412;
    hess[63] = t462;
    hess[64] = t473;
    hess[65] = t119
        * (-t107 - t109 + 2 * t110 * t180 * t27 * t94
           + 4 * t118 * t16 * t27 * t398 - t28 * t396 * t94 - t401 - t410);
    hess[66] = t486;
    hess[67] = t488;
    hess[68] = t490;
    hess[69] = t491;
    hess[70] = t492;
    hess[71] = t493;
    hess[72] = t233;
    hess[73] = t357;
    hess[74] = t429;
    hess[75] = t465;
    hess[76] = t476;
    hess[77] = t486;
    hess[78] = t119
        * (-t100 + t102 * t105 * t113 * t94 + t102 * t136 * t16 * t94 - t103
           - t119 * t17 * t231 - t124 * t499 - 4 * t149 * t499 - t494 - t496
           - t497 + t51 + t56 + t59 + t60 + t61 + t62 - t98);
    hess[79] = t506;
    hess[80] = t510;
    hess[81] = t512;
    hess[82] = t518;
    hess[83] = t522;
    hess[84] = t253;
    hess[85] = t366;
    hess[86] = t439;
    hess[87] = t467;
    hess[88] = t478;
    hess[89] = t488;
    hess[90] = t506;
    hess[91] = t119
        * (-t119 * t20 * t250 - t124 * t527 + t524 * t94 - t525 * t94
           + t527 * t528 + t531);
    hess[92] = t541;
    hess[93] = t542;
    hess[94] = t544;
    hess[95] = t547;
    hess[96] = t272;
    hess[97] = t380;
    hess[98] = t448;
    hess[99] = t469;
    hess[100] = t480;
    hess[101] = t490;
    hess[102] = t510;
    hess[103] = t541;
    hess[104] = t119
        * (-t119 * t23 * t269 - t124 * t552 + t528 * t552 + t549 * t94
           - t550 * t94 + t555);
    hess[105] = t556;
    hess[106] = t557;
    hess[107] = t558;
    hess[108] = t285;
    hess[109] = t385;
    hess[110] = t452;
    hess[111] = t470;
    hess[112] = t481;
    hess[113] = t491;
    hess[114] = t512;
    hess[115] = t542;
    hess[116] = t556;
    hess[117] = t560
        * (t102 * t105 * t113 + t102 * t136 * t16 - 2 * t216 * t274
           - t228 * t231 - 4 * t498 * t503 - t498 * t559);
    hess[118] = t562;
    hess[119] = t563;
    hess[120] = t291;
    hess[121] = t390;
    hess[122] = t457;
    hess[123] = t471;
    hess[124] = t482;
    hess[125] = t492;
    hess[126] = t518;
    hess[127] = t544;
    hess[128] = t557;
    hess[129] = t562;
    hess[130] = t560
        * (t105 * t113 * t523 + 2 * t105 * t21 * t239
           + 4 * t16 * t27 * t526 * t94 - t250 * t354 - t525 - t526 * t559);
    hess[131] = t564;
    hess[132] = t296;
    hess[133] = t395;
    hess[134] = t460;
    hess[135] = t472;
    hess[136] = t483;
    hess[137] = t493;
    hess[138] = t522;
    hess[139] = t547;
    hess[140] = t558;
    hess[141] = t563;
    hess[142] = t564;
    hess[143] = t560
        * (t105 * t113 * t548 + 2 * t105 * t24 * t260
           + 4 * t16 * t27 * t551 * t94 - t269 * t426 - t550 - t551 * t559);
}

// hess is (144×1) flattened in column-major order
void triangle_closest_point_hessian_0(
    double p_x,
    double p_y,
    double p_z,
    double t0_x,
    double t0_y,
    double t0_z,
    double t1_x,
    double t1_y,
    double t1_z,
    double t2_x,
    double t2_y,
    double t2_z,
    double hess[144])
{
    const auto t0 = t0_x * t1_x;
    const auto t1 = t0_y * t1_y;
    const auto t2 = 2 * t1;
    const auto t3 = t0_y * t2_y;
    const auto t4 = 2 * t3;
    const auto t5 = t0_x * t2_x;
    const auto t6 = 2 * t5;
    const auto t7 = t0_z * t1_z;
    const auto t8 = 2 * t7;
    const auto t9 = t0_z * t2_z;
    const auto t10 = 2 * t9;
    const auto t11 = t1_y * t2_y;
    const auto t12 = 2 * t11;
    const auto t13 = t1_z * t2_z;
    const auto t14 = 2 * t13;
    const auto t15 = t1_x * t2_x;
    const auto t16 = 2 * t15;
    const auto t17 = std::pow(t0_x, 2);
    const auto t18 = std::pow(t1_y, 2);
    const auto t19 = std::pow(t1_z, 2);
    const auto t20 = std::pow(t2_y, 2);
    const auto t21 = std::pow(t2_z, 2);
    const auto t22 = std::pow(t0_y, 2);
    const auto t23 = std::pow(t1_x, 2);
    const auto t24 = std::pow(t2_x, 2);
    const auto t25 = std::pow(t0_z, 2);
    const auto t26 = 2 * t0;
    const auto t27 = t0 * t10 + t0 * t12 + t0 * t14 - t0 * t2 + t0 * t4
        - t0 * t8 + t1 * t10 + t1 * t14 + t1 * t16 + t1 * t6 + t10 * t11
        + t10 * t15 - t10 * t18 - t10 * t23 + t11 * t6 - t12 * t13 - t12 * t15
        - t12 * t17 - t12 * t25 + t12 * t7 + t13 * t4 + t13 * t6 - t14 * t15
        - t14 * t17 - t14 * t22 + t15 * t4 - t16 * t22 - t16 * t25 + t16 * t7
        + t17 * t18 + t17 * t19 + t17 * t20 + t17 * t21 + t18 * t21 + t18 * t24
        + t18 * t25 - t18 * t6 + t19 * t20 + t19 * t22 + t19 * t24 - t19 * t4
        - t19 * t6 - t2 * t21 - t2 * t24 - t2 * t7 + t20 * t23 + t20 * t25
        - t20 * t26 - t20 * t8 + t21 * t22 + t21 * t23 - t21 * t26 + t22 * t23
        + t22 * t24 + t23 * t25 - t23 * t4 + t24 * t25 - t24 * t8 - t3 * t6
        + t4 * t7 - t4 * t9 + t6 * t7 - t6 * t9;
    const auto t28 = 1.0 / t27;
    const auto t29 = t0_x - t2_x;
    const auto t30 = 2 * t0_x;
    const auto t31 = -t30;
    const auto t32 = t2_x + t31;
    const auto t33 = t1_x + t32;
    const auto t34 = t0_x - t1_x;
    const auto t35 = t29 * t34;
    const auto t36 = t20 + t21;
    const auto t37 = -t10 + t25;
    const auto t38 = t22 - t4;
    const auto t39 = t36 + t37 + t38;
    const auto t40 = t17 - t6;
    const auto t41 = t24 + t39 + t40;
    const auto t42 = t0_x * t18;
    const auto t43 = t0_x * t19;
    const auto t44 = t0_x * t20;
    const auto t45 = t0_x * t21;
    const auto t46 = t1_x * t20;
    const auto t47 = t1_x * t21;
    const auto t48 = t18 * t2_x;
    const auto t49 = t19 * t2_x;
    const auto t50 = t11 * t1_x;
    const auto t51 = t13 * t1_x;
    const auto t52 = t11 * t2_x;
    const auto t53 = t13 * t2_x;
    const auto t54 = t1 * t1_x;
    const auto t55 = t2_x * t3;
    const auto t56 = t1_x * t7;
    const auto t57 = t2_x * t9;
    const auto t58 = t1 * t2_x;
    const auto t59 = t2_x * t7;
    const auto t60 = t58 + t59;
    const auto t61 = t1_x * t3;
    const auto t62 = t1_x * t9;
    const auto t63 = t61 + t62;
    const auto t64 = -t11 * t30 - t13 * t30 + t42 + t43 + t44 + t45 - t46 - t47
        - t48 - t49 + t50 + t51 + t52 + t53 - t54 - t55 - t56 - t57 + t60 + t63;
    const auto t65 = t11 + t13;
    const auto t66 = -t9;
    const auto t67 = t25 + t66 - t7;
    const auto t68 = -t3;
    const auto t69 = -t1 + t22 + t68;
    const auto t70 = t65 + t67 + t69;
    const auto t71 = -t5;
    const auto t72 = -t0 + t17 + t71;
    const auto t73 = t15 + t70 + t72;
    const auto t74 = t64 * t73;
    const auto t75 = 2 * t29;
    const auto t76 = t28 * t75;
    const auto t77 = -t13;
    const auto t78 = t66 + t7;
    const auto t79 = t77 + t78;
    const auto t80 = -t11;
    const auto t81 = t1 + t68;
    const auto t82 = t80 + t81;
    const auto t83 = t36 + t79 + t82;
    const auto t84 = -t15;
    const auto t85 = t0 + t71;
    const auto t86 = t84 + t85;
    const auto t87 = t24 + t83 + t86;
    const auto t88 = t28
        * (2 * t28 * t34 * t41 * t64 - t29 * t33 - 2 * t35 - t74 * t76 - t87);
    const auto t89 = t0_y - t2_y;
    const auto t90 = t34 * t89;
    const auto t91 = 2 * t90;
    const auto t92 = 2 * t0_y;
    const auto t93 = -t92;
    const auto t94 = t2_y + t93;
    const auto t95 = t1_y + t94;
    const auto t96 = t0_y * t23;
    const auto t97 = t0_y * t19;
    const auto t98 = t0_y * t24;
    const auto t99 = t0_y * t21;
    const auto t100 = t1_y * t24;
    const auto t101 = t1_y * t21;
    const auto t102 = t23 * t2_y;
    const auto t103 = t19 * t2_y;
    const auto t104 = t15 * t1_y;
    const auto t105 = t15 * t2_y;
    const auto t106 = t13 * t1_y;
    const auto t107 = t13 * t2_y;
    const auto t108 = t0 * t1_y;
    const auto t109 = t2_y * t5;
    const auto t110 = t1_y * t7;
    const auto t111 = t2_y * t9;
    const auto t112 = t0 * t2_y;
    const auto t113 = t2_y * t7;
    const auto t114 = t112 + t113;
    const auto t115 = t1_y * t5 + t1_y * t9;
    const auto t116 = -t100 - t101 - t102 - t103 + t104 + t105 + t106 + t107
        - t108 - t109 - t110 - t111 + t114 + t115 - t13 * t92 - t15 * t92 + t96
        + t97 + t98 + t99;
    const auto t117 = t73 * t76;
    const auto t118 =
        t28 * (-t116 * t117 + 2 * t116 * t28 * t34 * t41 - t29 * t95 - t91);
    const auto t119 = t0_z - t2_z;
    const auto t120 = t119 * t34;
    const auto t121 = 2 * t120;
    const auto t122 = 2 * t0_z;
    const auto t123 = -t122;
    const auto t124 = t123 + t2_z;
    const auto t125 = t124 + t1_z;
    const auto t126 = t0_z * t23;
    const auto t127 = t0_z * t18;
    const auto t128 = t0_z * t24;
    const auto t129 = t0_z * t20;
    const auto t130 = t1_z * t24;
    const auto t131 = t1_z * t20;
    const auto t132 = t23 * t2_z;
    const auto t133 = t18 * t2_z;
    const auto t134 = t15 * t1_z;
    const auto t135 = t15 * t2_z;
    const auto t136 = t11 * t1_z;
    const auto t137 = t11 * t2_z;
    const auto t138 = t0 * t1_z;
    const auto t139 = t2_z * t5;
    const auto t140 = t1 * t1_z;
    const auto t141 = t2_z * t3;
    const auto t142 = t0 * t2_z;
    const auto t143 = t1 * t2_z;
    const auto t144 = t142 + t143;
    const auto t145 = t1_z * t3 + t1_z * t5;
    const auto t146 = -t11 * t122 - t122 * t15 + t126 + t127 + t128 + t129
        - t130 - t131 - t132 - t133 + t134 + t135 + t136 + t137 - t138 - t139
        - t140 - t141 + t144 + t145;
    const auto t147 =
        t28 * (-t117 * t146 - t121 - t125 * t29 + 2 * t146 * t28 * t34 * t41);
    const auto t148 = t1_x * t22;
    const auto t149 = t1_x * t25;
    const auto t150 = t22 * t2_x;
    const auto t151 = t25 * t2_x;
    const auto t152 = t0_x * t3;
    const auto t153 = t0_x * t9;
    const auto t154 = t0_x * t1;
    const auto t155 = t0_x * t7;
    const auto t156 = t0_x * t11 + t0_x * t13;
    const auto t157 = t148 + t149 - t150 - t151 + t152 + t153 - t154 - t155
        + t156 - t44 - t45 + t46 + t47 - t52 - t53 + t55 + t57 + t60 - 2 * t61
        - 2 * t62;
    const auto t158 = t157 * t41;
    const auto t159 = 2 * t28;
    const auto t160 = t159 * t34;
    const auto t161 =
        t28 * (-t117 * t157 + t158 * t160 - std::pow(t29, 2) + t41);
    const auto t162 = t29 * t89;
    const auto t163 = t17 * t1_y;
    const auto t164 = t1_y * t25;
    const auto t165 = t17 * t2_y;
    const auto t166 = t25 * t2_y;
    const auto t167 = t0_y * t5;
    const auto t168 = t0_y * t9;
    const auto t169 = t0 * t0_y;
    const auto t170 = t0_y * t7;
    const auto t171 = t0_y * t13 + t0_y * t15;
    const auto t172 = -t10 * t1_y + t100 + t101 - t105 - t107 + t109 + t111
        + t114 + t163 + t164 - t165 - t166 + t167 + t168 - t169 - t170 + t171
        - t1_y * t6 - t98 - t99;
    const auto t173 = t28 * (-t117 * t172 - t162 + 2 * t172 * t28 * t34 * t41);
    const auto t174 = t119 * t29;
    const auto t175 = t17 * t1_z;
    const auto t176 = t1_z * t22;
    const auto t177 = t17 * t2_z;
    const auto t178 = t22 * t2_z;
    const auto t179 = t0_z * t5;
    const auto t180 = t0_z * t3;
    const auto t181 = t0 * t0_z;
    const auto t182 = t0_z * t1;
    const auto t183 = t0_z * t11 + t0_z * t15;
    const auto t184 = -t128 - t129 + t130 + t131 - t135 - t137 + t139 + t141
        + t144 + t175 + t176 - t177 - t178 + t179 + t180 - t181 - t182 + t183
        - t1_z * t4 - t1_z * t6;
    const auto t185 = t28 * (-t117 * t184 - t174 + 2 * t184 * t28 * t34 * t41);
    const auto t186 = -t148 - t149 + t150 + t151 - t152 - t153 + t154 + t155
        + t156 - t42 - t43 + t48 + t49 - t50 - t51 + t54 + t56 - 2 * t58
        - 2 * t59 + t63;
    const auto t187 = t160 * t41;
    const auto t188 =
        t0 + t1 - t17 - t22 - t25 + t3 + t5 + t7 + t77 + t80 + t84 + t9;
    const auto t189 = t28 * (-t117 * t186 + t186 * t187 + t188 + t35);
    const auto t190 = t0_y - t1_y;
    const auto t191 = t190 * t29;
    const auto t192 = t102 + t103 - t104 - t106 + t108 + t110 - 2 * t112
        - 2 * t113 + t115 - t163 - t164 + t165 + t166 - t167 - t168 + t169
        + t170 + t171 - t96 - t97;
    const auto t193 = t28 * (-t117 * t192 + t187 * t192 - t191 + t91);
    const auto t194 = t0_z - t1_z;
    const auto t195 = t194 * t29;
    const auto t196 = -t126 - t127 + t132 + t133 - t134 - t136 + t138 + t140
        - 2 * t142 - 2 * t143 + t145 - t175 - t176 + t177 + t178 - t179 - t180
        + t181 + t182 + t183;
    const auto t197 = t28 * (-t117 * t196 + t121 + t187 * t196 - t195);
    const auto t198 = 2 * t191;
    const auto t199 = t159 * t89;
    const auto t200 =
        t28 * (2 * t190 * t28 * t41 * t64 - t198 - t199 * t74 - t33 * t89);
    const auto t201 = t190 * t89;
    const auto t202 = t199 * t73;
    const auto t203 = t28
        * (2 * t116 * t190 * t28 * t41 - t116 * t202 - 2 * t201 - t87
           - t89 * t95);
    const auto t204 = t119 * t190;
    const auto t205 = 2 * t204;
    const auto t206 =
        t28 * (-t125 * t89 + 2 * t146 * t190 * t28 * t41 - t146 * t202 - t205);
    const auto t207 = t28 * (2 * t157 * t190 * t28 * t41 - t157 * t202 - t162);
    const auto t208 = t159 * t190;
    const auto t209 = t208 * t41;
    const auto t210 =
        t28 * (-t172 * t202 + t172 * t209 + t41 - std::pow(t89, 2));
    const auto t211 = t119 * t89;
    const auto t212 = t28 * (2 * t184 * t190 * t28 * t41 - t184 * t202 - t211);
    const auto t213 =
        t28 * (2 * t186 * t190 * t28 * t41 - t186 * t202 + t198 - t90);
    const auto t214 = t28 * (t188 - t192 * t202 + t192 * t209 + t201);
    const auto t215 = t194 * t89;
    const auto t216 = t28 * (-t196 * t202 + t196 * t209 + t205 - t215);
    const auto t217 = 2 * t195;
    const auto t218 = t119 * t159;
    const auto t219 =
        t28 * (-t119 * t33 + 2 * t194 * t28 * t41 * t64 - t217 - t218 * t74);
    const auto t220 = 2 * t215;
    const auto t221 = t218 * t73;
    const auto t222 =
        t28 * (2 * t116 * t194 * t28 * t41 - t116 * t221 - t119 * t95 - t220);
    const auto t223 = t119 * t194;
    const auto t224 = t28
        * (-t119 * t125 + 2 * t146 * t194 * t28 * t41 - t146 * t221 - 2 * t223
           - t87);
    const auto t225 = t28 * (2 * t157 * t194 * t28 * t41 - t157 * t221 - t174);
    const auto t226 = t28 * (2 * t172 * t194 * t28 * t41 - t172 * t221 - t211);
    const auto t227 = t159 * t194;
    const auto t228 = t227 * t41;
    const auto t229 =
        t28 * (-std::pow(t119, 2) - t184 * t221 + t184 * t228 + t41);
    const auto t230 =
        t28 * (-t120 + 2 * t186 * t194 * t28 * t41 - t186 * t221 + t217);
    const auto t231 =
        t28 * (2 * t192 * t194 * t28 * t41 - t192 * t221 - t204 + t220);
    const auto t232 = t28 * (t188 - t196 * t221 + t196 * t228 + t223);
    const auto t233 = p_x + t32;
    const auto t234 = p_x + t1_x + t31;
    const auto t235 = t234 * t75;
    const auto t236 = p_x - t0_x;
    const auto t237 = t236 * t34;
    const auto t238 = p_y - t0_y;
    const auto t239 = t190 * t238;
    const auto t240 = p_z - t0_z;
    const auto t241 = t194 * t240;
    const auto t242 = t239 + t241;
    const auto t243 = t237 + t242;
    const auto t244 = t243 * t41;
    const auto t245 = -t12;
    const auto t246 = -t14;
    const auto t247 = t18 + t19;
    const auto t248 = t28 * (t245 + t246 + t247 + t36);
    const auto t249 = t236 * t29;
    const auto t250 = t238 * t89;
    const auto t251 = t119 * t240;
    const auto t252 = t250 + t251;
    const auto t253 = t249 + t252;
    const auto t254 = t253 * t73;
    const auto t255 = std::pow(t27, -2);
    const auto t256 = 4 * t255;
    const auto t257 = t256 * std::pow(t64, 2);
    const auto t258 = t159 * t233;
    const auto t259 = t159 * t64;
    const auto t260 = t234 * t41;
    const auto t261 = t253 * t33;
    const auto t262 = 4 * t28;
    const auto t263 = t243 * t262;
    const auto t264 = t263 * t29;
    const auto t265 = t264 * t64;
    const auto t266 = -t237 - t239 - t241 + t253 + t87;
    const auto t267 = p_y + t94;
    const auto t268 = p_y + t1_y + t93;
    const auto t269 = t1_x * t1_y;
    const auto t270 = t2_x * t2_y;
    const auto t271 = t1_x * t2_y;
    const auto t272 = -t271;
    const auto t273 = t1_y * t2_x;
    const auto t274 = -t273;
    const auto t275 = t269 + t270 + t272 + t274;
    const auto t276 = t258 * t73;
    const auto t277 = t159 * t74;
    const auto t278 = 8 * t255;
    const auto t279 = t278 * t64;
    const auto t280 = 2 * t234;
    const auto t281 = t263 * t64;
    const auto t282 = t280 * t89 - t281 * t89;
    const auto t283 = -t116 * t264 + t268 * t75;
    const auto t284 = t28
        * (2 * t116 * t234 * t28 * t41 - t116 * t244 * t279
           + 8 * t116 * t253 * t255 * t64 * t73 + 2 * t116 * t253 * t28 * t33
           - t116 * t276 - t159 * t244 * t275 - t233 * t95
           + 2 * t253 * t275 * t28 * t73 + 2 * t253 * t28 * t64 * t95
           - t267 * t277 - t267 * t33 + 2 * t268 * t28 * t41 * t64 - t282
           - t283);
    const auto t285 = p_z + t124;
    const auto t286 = p_z + t123 + t1_z;
    const auto t287 = t1_x * t1_z;
    const auto t288 = t2_x * t2_z;
    const auto t289 = t1_x * t2_z;
    const auto t290 = -t289;
    const auto t291 = t1_z * t2_x;
    const auto t292 = -t291;
    const auto t293 = t287 + t288 + t290 + t292;
    const auto t294 = t119 * t280 - t119 * t281;
    const auto t295 = -t146 * t264 + t286 * t75;
    const auto t296 = t28
        * (-t125 * t233 + 2 * t125 * t253 * t28 * t64
           + 2 * t146 * t234 * t28 * t41 - t146 * t244 * t279
           + 8 * t146 * t253 * t255 * t64 * t73 + 2 * t146 * t253 * t28 * t33
           - t146 * t276 - t159 * t244 * t293 + 2 * t253 * t28 * t293 * t73
           - t277 * t285 + 2 * t28 * t286 * t41 * t64 - t285 * t33 - t294
           - t295);
    const auto t297 = -t249;
    const auto t298 = t236 * t41;
    const auto t299 = -t157 * t264 + t252;
    const auto t300 = t28
        * (2 * t157 * t234 * t28 * t41 - t157 * t244 * t279
           + 8 * t157 * t253 * t255 * t64 * t73 + 2 * t157 * t253 * t28 * t33
           - t157 * t276 - t159 * t244 * t83 - t233 * t29
           + 2 * t253 * t28 * t29 * t64 + 2 * t253 * t28 * t73 * t83
           - t259 * t298 - t297 - t299 - t41);
    const auto t301 = t238 * t29;
    const auto t302 = 2 * t301;
    const auto t303 = t0_y * t1_x;
    const auto t304 = -t303;
    const auto t305 = t1_y * t30;
    const auto t306 = 2 * t273;
    const auto t307 = t271 - t306;
    const auto t308 = t0_y * t2_x;
    const auto t309 = t2_y * t30;
    const auto t310 = t308 - t309;
    const auto t311 = t159 * (t270 + t304 + t305 + t307 + t310);
    const auto t312 = t238 * t41;
    const auto t313 = t159 * t260;
    const auto t314 = t172 * t264;
    const auto t315 = t253 * t64;
    const auto t316 = t172 * t253;
    const auto t317 = t159 * t33;
    const auto t318 = t172 * t279;
    const auto t319 = t28
        * (-t172 * t276 + t172 * t313 + t199 * t315 - t233 * t89 + t244 * t311
           - t244 * t318 - t254 * t311 + t254 * t318 - t259 * t312 + t302 + t314
           + t316 * t317);
    const auto t320 = t240 * t29;
    const auto t321 = 2 * t320;
    const auto t322 = t0_z * t1_x;
    const auto t323 = -t322;
    const auto t324 = t1_z * t30;
    const auto t325 = 2 * t291;
    const auto t326 = t289 - t325;
    const auto t327 = t0_z * t2_x;
    const auto t328 = t2_z * t30;
    const auto t329 = t327 - t328;
    const auto t330 = t159 * (t288 + t323 + t324 + t326 + t329);
    const auto t331 = t240 * t41;
    const auto t332 = t184 * t264;
    const auto t333 = t184 * t253;
    const auto t334 = t184 * t279;
    const auto t335 = t28
        * (-t119 * t233 - t184 * t276 + t184 * t313 + t218 * t315 + t244 * t330
           - t244 * t334 - t254 * t330 + t254 * t334 - t259 * t331 + t317 * t333
           + t321 + t332);
    const auto t336 = t186 * t264;
    const auto t337 = -t19 + t78;
    const auto t338 = -t18 + t81;
    const auto t339 = t159 * (t337 + t338 + t65);
    const auto t340 = t159 * t261;
    const auto t341 = t186 * t279;
    const auto t342 = -t250;
    const auto t343 = -t251;
    const auto t344 = 2 * t237 + 2 * t239 + 2 * t241 + t297 + t342 + t343 + t73;
    const auto t345 = t28
        * (t160 * t315 - t186 * t276 + t186 * t313 + t186 * t340 - t233 * t34
           + t235 + t236 * t277 + t236 * t33 + t244 * t339 - t244 * t341
           - t254 * t339 + t254 * t341 - t265 + t336 + t344);
    const auto t346 = -t308;
    const auto t347 = 2 * t271;
    const auto t348 = t273 - t347;
    const auto t349 = t303 - t305;
    const auto t350 = t159 * (t269 + t309 + t346 + t348 + t349);
    const auto t351 = t192 * t264;
    const auto t352 = t192 * t279;
    const auto t353 = t28
        * (-t190 * t233 - t192 * t276 + t192 * t313 + t192 * t340 + t208 * t315
           + t238 * t277 + t238 * t33 + t244 * t350 - t244 * t352 - t254 * t350
           + t254 * t352 + t282 + t351);
    const auto t354 = -t327;
    const auto t355 = 2 * t289;
    const auto t356 = t291 - t355;
    const auto t357 = t322 - t324;
    const auto t358 = t159 * (t287 + t328 + t354 + t356 + t357);
    const auto t359 = t196 * t264;
    const auto t360 = t196 * t279;
    const auto t361 = t28
        * (-t194 * t233 - t196 * t276 + t196 * t313 + t196 * t340 + t227 * t315
           + t240 * t277 + t240 * t33 + t244 * t358 - t244 * t360 - t254 * t358
           + t254 * t360 + t294 + t359);
    const auto t362 = 2 * t89;
    const auto t363 = t268 * t362;
    const auto t364 = -t16;
    const auto t365 = t21 + t24;
    const auto t366 = t19 + t23;
    const auto t367 = t28 * (t246 + t364 + t365 + t366);
    const auto t368 = std::pow(t116, 2) * t256;
    const auto t369 = t116 * t159;
    const auto t370 = t267 * t73;
    const auto t371 = t268 * t41;
    const auto t372 = t253 * t95;
    const auto t373 = t263 * t89;
    const auto t374 = t116 * t373;
    const auto t375 = t1_y * t1_z;
    const auto t376 = t2_y * t2_z;
    const auto t377 = t1_y * t2_z;
    const auto t378 = -t377;
    const auto t379 = t1_z * t2_y;
    const auto t380 = -t379;
    const auto t381 = t375 + t376 + t378 + t380;
    const auto t382 = t146 * t159;
    const auto t383 = t369 * t73;
    const auto t384 = t116 * t278;
    const auto t385 = 2 * t119;
    const auto t386 = t119 * t263;
    const auto t387 = -t116 * t386 + t268 * t385;
    const auto t388 = -t146 * t373 + t286 * t362;
    const auto t389 = t28
        * (2 * t116 * t125 * t253 * t28 + 8 * t116 * t146 * t253 * t255 * t73
           + 2 * t116 * t28 * t286 * t41 - t125 * t267 - t146 * t244 * t384
           + 2 * t146 * t253 * t28 * t95 + 2 * t146 * t268 * t28 * t41
           - t159 * t244 * t381 + 2 * t253 * t28 * t381 * t73 - t285 * t383
           - t285 * t95 - t370 * t382 - t387 - t388);
    const auto t390 = t236 * t89;
    const auto t391 = 2 * t390;
    const auto t392 = t0_x * t1_y;
    const auto t393 = -t392;
    const auto t394 = t1_x * t92;
    const auto t395 = t0_x * t2_y;
    const auto t396 = t2_x * t92;
    const auto t397 = t395 - t396;
    const auto t398 = t159 * (t270 + t348 + t393 + t394 + t397);
    const auto t399 = t158 * t159;
    const auto t400 = t116 * t253;
    const auto t401 = t157 * t373;
    const auto t402 = t157 * t253;
    const auto t403 = t159 * t95;
    const auto t404 = t159 * t370;
    const auto t405 = t157 * t384;
    const auto t406 = t28
        * (-t157 * t404 + t244 * t398 - t244 * t405 - t254 * t398 + t254 * t405
           - t267 * t29 + t268 * t399 - t298 * t369 + t391 + t400 * t76 + t401
           + t402 * t403);
    const auto t407 = t365 + t79 + t86;
    const auto t408 = -t172 * t373 + t249 + t251;
    const auto t409 = t28
        * (8 * t116 * t172 * t253 * t255 * t73 + 2 * t116 * t253 * t28 * t89
           - t159 * t244 * t407 - t172 * t244 * t384
           + 2 * t172 * t253 * t28 * t95 + 2 * t172 * t268 * t28 * t41
           - t172 * t404 + 2 * t253 * t28 * t407 * t73 - t267 * t89
           - t312 * t369 - t342 - t408 - t41);
    const auto t410 = t240 * t89;
    const auto t411 = 2 * t410;
    const auto t412 = t0_z * t1_y;
    const auto t413 = -t412;
    const auto t414 = t1_z * t92;
    const auto t415 = 2 * t379;
    const auto t416 = t377 - t415;
    const auto t417 = t0_z * t2_y;
    const auto t418 = t2_z * t92;
    const auto t419 = t417 - t418;
    const auto t420 = t159 * (t376 + t413 + t414 + t416 + t419);
    const auto t421 = t159 * t371;
    const auto t422 = t184 * t373;
    const auto t423 = t184 * t384;
    const auto t424 = t28
        * (-t119 * t267 - t184 * t404 + t184 * t421 + t218 * t400 + t244 * t420
           - t244 * t423 - t254 * t420 + t254 * t423 - t331 * t369 + t333 * t403
           + t411 + t422);
    const auto t425 = -t395;
    const auto t426 = t392 - t394;
    const auto t427 = t159 * (t269 + t307 + t396 + t425 + t426);
    const auto t428 = t159 * t372;
    const auto t429 = t186 * t373;
    const auto t430 = t186 * t384;
    const auto t431 = t28
        * (t160 * t400 - t186 * t404 + t186 * t421 + t186 * t428 + t236 * t383
           + t236 * t95 + t244 * t427 - t244 * t430 - t254 * t427 + t254 * t430
           - t267 * t34 + t283 + t429);
    const auto t432 = t192 * t373;
    const auto t433 = t13 + t15;
    const auto t434 = -t23 + t85;
    const auto t435 = t159 * (t337 + t433 + t434);
    const auto t436 = t192 * t384;
    const auto t437 = t28
        * (-t190 * t267 - t192 * t404 + t192 * t421 + t192 * t428 + t208 * t400
           + t238 * t383 + t238 * t95 + t244 * t435 - t244 * t436 - t254 * t435
           + t254 * t436 + t344 + t363 - t374 + t432);
    const auto t438 = -t417;
    const auto t439 = 2 * t377;
    const auto t440 = t379 - t439;
    const auto t441 = t412 - t414;
    const auto t442 = t159 * (t375 + t418 + t438 + t440 + t441);
    const auto t443 = t196 * t373;
    const auto t444 = t196 * t384;
    const auto t445 = t28
        * (-t194 * t267 - t196 * t404 + t196 * t421 + t196 * t428 + t227 * t400
           + t240 * t383 + t240 * t95 + t244 * t442 - t244 * t444 - t254 * t442
           + t254 * t444 + t387 + t443);
    const auto t446 = t286 * t385;
    const auto t447 = t20 + t24;
    const auto t448 = t18 + t23;
    const auto t449 = t28 * (t245 + t364 + t447 + t448);
    const auto t450 = std::pow(t146, 2) * t256;
    const auto t451 = t285 * t73;
    const auto t452 = t286 * t41;
    const auto t453 = t125 * t253;
    const auto t454 = t146 * t386;
    const auto t455 = t119 * t236;
    const auto t456 = 2 * t455;
    const auto t457 = t0_x * t1_z;
    const auto t458 = -t457;
    const auto t459 = t122 * t1_x;
    const auto t460 = t0_x * t2_z;
    const auto t461 = t122 * t2_x;
    const auto t462 = t460 - t461;
    const auto t463 = t159 * (t288 + t356 + t458 + t459 + t462);
    const auto t464 = t146 * t253;
    const auto t465 = t157 * t386;
    const auto t466 = t125 * t159;
    const auto t467 = t159 * t451;
    const auto t468 = t146 * t278;
    const auto t469 = t157 * t244;
    const auto t470 = t254 * t468;
    const auto t471 = t28
        * (-t157 * t467 + t157 * t470 + t244 * t463 - t254 * t463 - t285 * t29
           + t286 * t399 - t298 * t382 + t402 * t466 + t456 + t464 * t76 + t465
           - t468 * t469);
    const auto t472 = t119 * t238;
    const auto t473 = 2 * t472;
    const auto t474 = t0_y * t1_z;
    const auto t475 = -t474;
    const auto t476 = t122 * t1_y;
    const auto t477 = t0_y * t2_z;
    const auto t478 = t122 * t2_y;
    const auto t479 = t477 - t478;
    const auto t480 = t159 * (t376 + t440 + t475 + t476 + t479);
    const auto t481 = t159 * t452;
    const auto t482 = t172 * t386;
    const auto t483 = t244 * t468;
    const auto t484 = t28
        * (-t172 * t467 + t172 * t470 + t172 * t481 - t172 * t483 + t199 * t464
           + t244 * t480 - t254 * t480 - t285 * t89 - t312 * t382 + t316 * t466
           + t473 + t482);
    const auto t485 = t447 + t82 + t86;
    const auto t486 = -t184 * t386 + t249 + t250;
    const auto t487 = t28
        * (2 * t119 * t146 * t253 * t28 - t119 * t285
           + 2 * t125 * t184 * t253 * t28 + 8 * t146 * t184 * t253 * t255 * t73
           - t159 * t244 * t485 + 2 * t184 * t28 * t286 * t41 - t184 * t467
           - t184 * t483 + 2 * t253 * t28 * t485 * t73 - t331 * t382 - t343
           - t41 - t486);
    const auto t488 = -t460;
    const auto t489 = t457 - t459;
    const auto t490 = t159 * (t287 + t326 + t461 + t488 + t489);
    const auto t491 = t382 * t73;
    const auto t492 = t159 * t453;
    const auto t493 = t186 * t386;
    const auto t494 = t28
        * (t125 * t236 + t160 * t464 - t186 * t467 + t186 * t470 + t186 * t481
           - t186 * t483 + t186 * t492 + t236 * t491 + t244 * t490 - t254 * t490
           - t285 * t34 + t295 + t493);
    const auto t495 = -t477;
    const auto t496 = t474 - t476;
    const auto t497 = t159 * (t375 + t416 + t478 + t495 + t496);
    const auto t498 = t192 * t386;
    const auto t499 = t28
        * (t125 * t238 - t190 * t285 - t192 * t467 + t192 * t470 + t192 * t481
           - t192 * t483 + t192 * t492 + t208 * t464 + t238 * t491 + t244 * t497
           - t254 * t497 + t388 + t498);
    const auto t500 = t196 * t386;
    const auto t501 = t11 + t15;
    const auto t502 = t159 * (t338 + t434 + t501);
    const auto t503 = t28
        * (t125 * t240 - t194 * t285 - t196 * t467 + t196 * t470 + t196 * t481
           - t196 * t483 + t196 * t492 + t227 * t464 + t240 * t491 + t244 * t502
           - t254 * t502 + t344 + t446 - t454 + t500);
    const auto t504 = std::pow(t157, 2);
    const auto t505 = 2 * t255;
    const auto t506 = t0_x * t0_y;
    const auto t507 = t270 + t346 + t425 + t506;
    const auto t508 = t172 * t262;
    const auto t509 = t505
        * (4 * t157 * t172 * t253 * t28 * t73 + t157 * t253 * t89 - t158 * t238
           + t172 * t253 * t29 - t172 * t298 - t244 * t507 + t253 * t507 * t73
           - t469 * t508);
    const auto t510 = t0_x * t0_z;
    const auto t511 = t288 + t354 + t488 + t510;
    const auto t512 = t505
        * (t119 * t157 * t253 + 4 * t157 * t184 * t253 * t28 * t73 - t158 * t240
           + t184 * t253 * t29 - t184 * t262 * t469 - t184 * t298 - t244 * t511
           + t253 * t511 * t73);
    const auto t513 = t159 * t298;
    const auto t514 = t159 * t70;
    const auto t515 = t159 * t73;
    const auto t516 = t157 * t515;
    const auto t517 = t253 * t76;
    const auto t518 = t186 * t278;
    const auto t519 = t157 * t254;
    const auto t520 = t28
        * (t160 * t402 - t186 * t513 + t186 * t517 + t236 * t516 - t244 * t514
           + t254 * t514 + t299 - t469 * t518 + t518 * t519);
    const auto t521 = t159 * (t274 + t310 + t347 + t426 + t506);
    const auto t522 = t192 * t278;
    const auto t523 = t28
        * (-t192 * t513 + t192 * t517 + t208 * t402 + t238 * t516 + t244 * t521
           - t254 * t521 + t301 - t391 - t401 - t469 * t522 + t519 * t522);
    const auto t524 = t159 * (t292 + t329 + t355 + t489 + t510);
    const auto t525 = t196 * t278;
    const auto t526 = t28
        * (-t196 * t513 + t196 * t517 + t227 * t402 + t240 * t516 + t244 * t524
           - t254 * t524 + t320 - t456 - t465 - t469 * t525 + t519 * t525);
    const auto t527 = t365 + t37 + t40;
    const auto t528 = std::pow(t172, 2);
    const auto t529 = t0_y * t0_z;
    const auto t530 = t376 + t438 + t495 + t529;
    const auto t531 = t505
        * (t119 * t172 * t253 + 4 * t172 * t184 * t253 * t28 * t73 - t172 * t331
           - t184 * t244 * t508 + t184 * t253 * t89 - t184 * t312 - t244 * t530
           + t253 * t530 * t73);
    const auto t532 = t159 * (t272 + t306 + t349 + t397 + t506);
    const auto t533 = t159 * t312;
    const auto t534 = t199 * t253;
    const auto t535 = t172 * t515;
    const auto t536 = t172 * t518;
    const auto t537 = t28
        * (t160 * t316 - t186 * t533 + t186 * t534 + t236 * t535 + t244 * t532
           - t244 * t536 - t254 * t532 + t254 * t536 - t302 - t314 + t390);
    const auto t538 = t159 * (t433 + t67 + t72);
    const auto t539 = t172 * t522;
    const auto t540 = t28
        * (-t192 * t533 + t192 * t534 + t208 * t316 + t238 * t535 - t244 * t538
           - t244 * t539 + t254 * t538 + t254 * t539 + t408);
    const auto t541 = t159 * (t380 + t419 + t439 + t496 + t529);
    const auto t542 = t172 * t525;
    const auto t543 = t28
        * (-t196 * t533 + t196 * t534 + t227 * t316 + t240 * t535 + t244 * t541
           - t244 * t542 - t254 * t541 + t254 * t542 + t410 - t473 - t482);
    const auto t544 = t38 + t40 + t447;
    const auto t545 = std::pow(t184, 2);
    const auto t546 = t159 * (t290 + t325 + t357 + t462 + t510);
    const auto t547 = t159 * t331;
    const auto t548 = t218 * t253;
    const auto t549 = t184 * t515;
    const auto t550 = t184 * t518;
    const auto t551 = t28
        * (t160 * t333 - t186 * t547 + t186 * t548 + t236 * t549 + t244 * t546
           - t244 * t550 - t254 * t546 + t254 * t550 - t321 - t332 + t455);
    const auto t552 = t159 * (t378 + t415 + t441 + t479 + t529);
    const auto t553 = t184 * t522;
    const auto t554 = t28
        * (-t192 * t547 + t192 * t548 + t208 * t333 + t238 * t549 + t244 * t552
           - t244 * t553 - t254 * t552 + t254 * t553 - t411 - t422 + t472);
    const auto t555 = t159 * (t501 + t69 + t72);
    const auto t556 = t184 * t525;
    const auto t557 = t28
        * (-t196 * t547 + t196 * t548 + t227 * t333 + t240 * t549 - t244 * t555
           - t244 * t556 + t254 * t555 + t254 * t556 + t486);
    const auto t558 = t25 - t8;
    const auto t559 = -t2 + t22;
    const auto t560 = t247 + t558 + t559;
    const auto t561 = std::pow(t186, 2);
    const auto t562 = t159 * (t269 + t304 + t393 + t506);
    const auto t563 = t160 * t253;
    const auto t564 = t186 * t253;
    const auto t565 = t236 * t515;
    const auto t566 = t186 * t515;
    const auto t567 = t192 * t518;
    const auto t568 = t28
        * (t190 * t236 + t192 * t563 + t192 * t565 + t208 * t564 + t238 * t34
           + t238 * t566 - t244 * t562 - t244 * t567 + t254 * t562 + t254 * t567
           - t351 - t429);
    const auto t569 = t159 * (t287 + t323 + t458 + t510);
    const auto t570 = t196 * t518;
    const auto t571 = t28
        * (t194 * t236 + t196 * t563 + t196 * t565 + t227 * t564 + t240 * t34
           + t240 * t566 - t244 * t569 - t244 * t570 + t254 * t569 + t254 * t570
           - t359 - t493);
    const auto t572 = t17 - t26;
    const auto t573 = t366 + t558 + t572;
    const auto t574 = std::pow(t192, 2);
    const auto t575 = t159 * (t375 + t413 + t475 + t529);
    const auto t576 = t196 * t522;
    const auto t577 = t28
        * (t190 * t240 + t192 * t227 * t253 + t192 * t240 * t515 + t194 * t238
           + t196 * t208 * t253 + t196 * t238 * t515 - t244 * t575 - t244 * t576
           + t254 * t575 + t254 * t576 - t443 - t498);
    const auto t578 = t448 + t559 + t572;
    const auto t579 = std::pow(t196, 2);
    hess[0] = 0;
    hess[1] = 0;
    hess[2] = 0;
    hess[3] = t88;
    hess[4] = t118;
    hess[5] = t147;
    hess[6] = t161;
    hess[7] = t173;
    hess[8] = t185;
    hess[9] = t189;
    hess[10] = t193;
    hess[11] = t197;
    hess[12] = 0;
    hess[13] = 0;
    hess[14] = 0;
    hess[15] = t200;
    hess[16] = t203;
    hess[17] = t206;
    hess[18] = t207;
    hess[19] = t210;
    hess[20] = t212;
    hess[21] = t213;
    hess[22] = t214;
    hess[23] = t216;
    hess[24] = 0;
    hess[25] = 0;
    hess[26] = 0;
    hess[27] = t219;
    hess[28] = t222;
    hess[29] = t224;
    hess[30] = t225;
    hess[31] = t226;
    hess[32] = t229;
    hess[33] = t230;
    hess[34] = t231;
    hess[35] = t232;
    hess[36] = t88;
    hess[37] = t200;
    hess[38] = t219;
    hess[39] = t159
        * (-t233 * t33 - t235 + t244 * t248 - t244 * t257 - t248 * t254
           + t254 * t257 - t258 * t74 + t259 * t260 + t259 * t261 + t265
           + t266);
    hess[40] = t284;
    hess[41] = t296;
    hess[42] = t300;
    hess[43] = t319;
    hess[44] = t335;
    hess[45] = t345;
    hess[46] = t353;
    hess[47] = t361;
    hess[48] = t118;
    hess[49] = t203;
    hess[50] = t222;
    hess[51] = t284;
    hess[52] = t159
        * (t244 * t367 - t244 * t368 - t254 * t367 + t254 * t368 + t266
           - t267 * t95 - t363 - t369 * t370 + t369 * t371 + t369 * t372
           + t374);
    hess[53] = t389;
    hess[54] = t406;
    hess[55] = t409;
    hess[56] = t424;
    hess[57] = t431;
    hess[58] = t437;
    hess[59] = t445;
    hess[60] = t147;
    hess[61] = t206;
    hess[62] = t224;
    hess[63] = t296;
    hess[64] = t389;
    hess[65] = t159
        * (-t125 * t285 + t244 * t449 - t244 * t450 - t254 * t449 + t254 * t450
           + t266 - t382 * t451 + t382 * t452 + t382 * t453 - t446 + t454);
    hess[66] = t471;
    hess[67] = t484;
    hess[68] = t487;
    hess[69] = t494;
    hess[70] = t499;
    hess[71] = t503;
    hess[72] = t161;
    hess[73] = t207;
    hess[74] = t225;
    hess[75] = t300;
    hess[76] = t406;
    hess[77] = t471;
    hess[78] = t505
        * (2 * t157 * t253 * t29 - 2 * t157 * t298 + t243 * t39 * t41
           - t244 * t262 * t504 + 4 * t253 * t28 * t504 * t73 - t254 * t39);
    hess[79] = t509;
    hess[80] = t512;
    hess[81] = t520;
    hess[82] = t523;
    hess[83] = t526;
    hess[84] = t173;
    hess[85] = t210;
    hess[86] = t226;
    hess[87] = t319;
    hess[88] = t409;
    hess[89] = t484;
    hess[90] = t509;
    hess[91] = t505
        * (2 * t172 * t253 * t89 - 2 * t172 * t312 + t243 * t41 * t527
           - t244 * t262 * t528 + 4 * t253 * t28 * t528 * t73 - t254 * t527);
    hess[92] = t531;
    hess[93] = t537;
    hess[94] = t540;
    hess[95] = t543;
    hess[96] = t185;
    hess[97] = t212;
    hess[98] = t229;
    hess[99] = t335;
    hess[100] = t424;
    hess[101] = t487;
    hess[102] = t512;
    hess[103] = t531;
    hess[104] = t505
        * (2 * t119 * t184 * t253 - 2 * t184 * t331 + t243 * t41 * t544
           - t244 * t262 * t545 + 4 * t253 * t28 * t545 * t73 - t254 * t544);
    hess[105] = t551;
    hess[106] = t554;
    hess[107] = t557;
    hess[108] = t189;
    hess[109] = t213;
    hess[110] = t230;
    hess[111] = t345;
    hess[112] = t431;
    hess[113] = t494;
    hess[114] = t520;
    hess[115] = t537;
    hess[116] = t551;
    hess[117] = t159
        * (2 * t186 * t236 * t28 * t73 + 2 * t186 * t253 * t28 * t34 - t242
           + t243 * t28 * t41 * t560 - t244 * t256 * t561
           + 4 * t253 * t255 * t561 * t73 - t254 * t28 * t560 - t336);
    hess[118] = t568;
    hess[119] = t571;
    hess[120] = t193;
    hess[121] = t214;
    hess[122] = t231;
    hess[123] = t353;
    hess[124] = t437;
    hess[125] = t499;
    hess[126] = t523;
    hess[127] = t540;
    hess[128] = t554;
    hess[129] = t568;
    hess[130] = t159
        * (2 * t190 * t192 * t253 * t28 + 2 * t192 * t238 * t28 * t73 - t237
           - t241 + t243 * t28 * t41 * t573 - t244 * t256 * t574
           + 4 * t253 * t255 * t574 * t73 - t254 * t28 * t573 - t432);
    hess[131] = t577;
    hess[132] = t197;
    hess[133] = t216;
    hess[134] = t232;
    hess[135] = t361;
    hess[136] = t445;
    hess[137] = t503;
    hess[138] = t526;
    hess[139] = t543;
    hess[140] = t557;
    hess[141] = t571;
    hess[142] = t577;
    hess[143] = t159
        * (2 * t194 * t196 * t253 * t28 + 2 * t196 * t240 * t28 * t73 - t237
           - t239 + t243 * t28 * t41 * t578 - t244 * t256 * t579
           + 4 * t253 * t255 * t579 * t73 - t254 * t28 * t578 - t500);
}

// hess is (144×1) flattened in column-major order
void triangle_closest_point_hessian_1(
    double p_x,
    double p_y,
    double p_z,
    double t0_x,
    double t0_y,
    double t0_z,
    double t1_x,
    double t1_y,
    double t1_z,
    double t2_x,
    double t2_y,
    double t2_z,
    double hess[144])
{
    const auto t0 = t0_y * t1_y;
    const auto t1 = t0_x * t1_x;
    const auto t2 = 2 * t1;
    const auto t3 = t0_y * t2_y;
    const auto t4 = t0_x * t2_x;
    const auto t5 = 2 * t0;
    const auto t6 = 2 * t3;
    const auto t7 = t0_z * t1_z;
    const auto t8 = t0_z * t2_z;
    const auto t9 = 2 * t7;
    const auto t10 = 2 * t8;
    const auto t11 = t1_y * t2_y;
    const auto t12 = t1_z * t2_z;
    const auto t13 = 2 * t11;
    const auto t14 = 2 * t12;
    const auto t15 = t1_x * t2_x;
    const auto t16 = 2 * t15;
    const auto t17 = std::pow(t0_x, 2);
    const auto t18 = std::pow(t1_y, 2);
    const auto t19 = std::pow(t1_z, 2);
    const auto t20 = std::pow(t2_y, 2);
    const auto t21 = std::pow(t2_z, 2);
    const auto t22 = std::pow(t0_y, 2);
    const auto t23 = std::pow(t1_x, 2);
    const auto t24 = std::pow(t2_x, 2);
    const auto t25 = std::pow(t0_z, 2);
    const auto t26 = 2 * t4;
    const auto t27 = -t0 * t2 - t10 * t18 - t10 * t23 - t10 * t4 + t11 * t2
        + t11 * t9 - t12 * t13 + t12 * t2 + t12 * t5 - t13 * t15 - t13 * t17
        - t13 * t25 + t13 * t4 + t13 * t8 - t14 * t15 - t14 * t17 - t14 * t22
        + t14 * t3 + t14 * t4 + t15 * t5 + t15 * t9 - t16 * t22 - t16 * t25
        + t16 * t3 + t16 * t8 + t17 * t18 + t17 * t19 + t17 * t20 + t17 * t21
        + t18 * t21 + t18 * t24 + t18 * t25 - t18 * t26 + t19 * t20 + t19 * t22
        + t19 * t24 - t19 * t26 - t19 * t6 - t2 * t20 - t2 * t21 + t2 * t3
        - t2 * t7 + t2 * t8 + t20 * t23 + t20 * t25 - t20 * t9 + t21 * t22
        + t21 * t23 - t21 * t5 + t22 * t23 + t22 * t24 + t23 * t25 - t23 * t6
        + t24 * t25 - t24 * t5 - t24 * t9 + t3 * t9 + t4 * t5 - t4 * t6
        + t4 * t9 - t5 * t7 + t5 * t8 - t6 * t8;
    const auto t28 = 1.0 / t27;
    const auto t29 = t0_x - t1_x;
    const auto t30 = 2 * t0_x;
    const auto t31 = -t30;
    const auto t32 = t1_x + t31;
    const auto t33 = t2_x + t32;
    const auto t34 = t0_x - t2_x;
    const auto t35 = t29 * t34;
    const auto t36 = t18 + t19;
    const auto t37 = t25 - t9;
    const auto t38 = t22 - t5;
    const auto t39 = t36 + t37 + t38;
    const auto t40 = t17 - t2;
    const auto t41 = t23 + t39 + t40;
    const auto t42 = t0_x * t18;
    const auto t43 = t0_x * t19;
    const auto t44 = t0_x * t20;
    const auto t45 = t0_x * t21;
    const auto t46 = t1_x * t20;
    const auto t47 = t1_x * t21;
    const auto t48 = t18 * t2_x;
    const auto t49 = t19 * t2_x;
    const auto t50 = t11 * t1_x;
    const auto t51 = t12 * t1_x;
    const auto t52 = t11 * t2_x;
    const auto t53 = t12 * t2_x;
    const auto t54 = t0 * t1_x;
    const auto t55 = t2_x * t3;
    const auto t56 = t1_x * t7;
    const auto t57 = t2_x * t8;
    const auto t58 = t0 * t2_x;
    const auto t59 = t2_x * t7;
    const auto t60 = t58 + t59;
    const auto t61 = t1_x * t3;
    const auto t62 = t1_x * t8;
    const auto t63 = t61 + t62;
    const auto t64 = -t11 * t30 - t12 * t30 + t42 + t43 + t44 + t45 - t46 - t47
        - t48 - t49 + t50 + t51 + t52 + t53 - t54 - t55 - t56 - t57 + t60 + t63;
    const auto t65 = -t0;
    const auto t66 = -t7;
    const auto t67 = t65 + t66;
    const auto t68 = -t8;
    const auto t69 = t12 + t25 + t68;
    const auto t70 = -t3;
    const auto t71 = t11 + t22 + t70;
    const auto t72 = t67 + t69 + t71;
    const auto t73 = -t4;
    const auto t74 = -t1;
    const auto t75 = t15 + t17 + t73 + t74;
    const auto t76 = t72 + t75;
    const auto t77 = t64 * t76;
    const auto t78 = 2 * t29;
    const auto t79 = t28 * t78;
    const auto t80 = -t15;
    const auto t81 = -t11;
    const auto t82 = -t12;
    const auto t83 = t81 + t82;
    const auto t84 = t3 + t4 + t8 + t80 + t83;
    const auto t85 = t23 + t36 + t67 + t74 + t84;
    const auto t86 = t28
        * (2 * t28 * t34 * t41 * t64 - t29 * t33 - 2 * t35 - t77 * t79 - t85);
    const auto t87 = 2 * t0_y;
    const auto t88 = -t87;
    const auto t89 = t1_y + t88;
    const auto t90 = t2_y + t89;
    const auto t91 = t0_y - t1_y;
    const auto t92 = t34 * t91;
    const auto t93 = 2 * t92;
    const auto t94 = t0_y * t23;
    const auto t95 = t0_y * t19;
    const auto t96 = t0_y * t24;
    const auto t97 = t0_y * t21;
    const auto t98 = t1_y * t24;
    const auto t99 = t1_y * t21;
    const auto t100 = t23 * t2_y;
    const auto t101 = t19 * t2_y;
    const auto t102 = t15 * t1_y;
    const auto t103 = t15 * t2_y;
    const auto t104 = t12 * t1_y;
    const auto t105 = t12 * t2_y;
    const auto t106 = t1 * t1_y;
    const auto t107 = t2_y * t4;
    const auto t108 = t1_y * t7;
    const auto t109 = t2_y * t8;
    const auto t110 = t1 * t2_y + t2_y * t7;
    const auto t111 = t1_y * t4;
    const auto t112 = t1_y * t8;
    const auto t113 = t111 + t112;
    const auto t114 = -t100 - t101 + t102 + t103 + t104 + t105 - t106 - t107
        - t108 - t109 + t110 + t113 - t12 * t87 - t15 * t87 + t94 + t95 + t96
        + t97 - t98 - t99;
    const auto t115 = t76 * t79;
    const auto t116 =
        t28 * (-t114 * t115 + 2 * t114 * t28 * t34 * t41 - t29 * t90 - t93);
    const auto t117 = 2 * t0_z;
    const auto t118 = -t117;
    const auto t119 = t118 + t1_z;
    const auto t120 = t119 + t2_z;
    const auto t121 = t0_z - t1_z;
    const auto t122 = t121 * t34;
    const auto t123 = 2 * t122;
    const auto t124 = t0_z * t23;
    const auto t125 = t0_z * t18;
    const auto t126 = t0_z * t24;
    const auto t127 = t0_z * t20;
    const auto t128 = t1_z * t24;
    const auto t129 = t1_z * t20;
    const auto t130 = t23 * t2_z;
    const auto t131 = t18 * t2_z;
    const auto t132 = t15 * t1_z;
    const auto t133 = t15 * t2_z;
    const auto t134 = t11 * t1_z;
    const auto t135 = t11 * t2_z;
    const auto t136 = t1 * t1_z;
    const auto t137 = t2_z * t4;
    const auto t138 = t0 * t1_z;
    const auto t139 = t2_z * t3;
    const auto t140 = t0 * t2_z + t1 * t2_z;
    const auto t141 = t1_z * t4;
    const auto t142 = t1_z * t3;
    const auto t143 = t141 + t142;
    const auto t144 = -t11 * t117 - t117 * t15 + t124 + t125 + t126 + t127
        - t128 - t129 - t130 - t131 + t132 + t133 + t134 + t135 - t136 - t137
        - t138 - t139 + t140 + t143;
    const auto t145 =
        t28 * (-t115 * t144 - t120 * t29 - t123 + 2 * t144 * t28 * t34 * t41);
    const auto t146 = t1_x * t22;
    const auto t147 = t1_x * t25;
    const auto t148 = t22 * t2_x;
    const auto t149 = t25 * t2_x;
    const auto t150 = t0_x * t3;
    const auto t151 = t0_x * t8;
    const auto t152 = t0 * t0_x;
    const auto t153 = t0_x * t7;
    const auto t154 = t0_x * t11 + t0_x * t12;
    const auto t155 = t146 + t147 - t148 - t149 + t150 + t151 - t152 - t153
        + t154 - t44 - t45 + t46 + t47 - t52 - t53 + t55 + t57 + t60 - 2 * t61
        - 2 * t62;
    const auto t156 = 2 * t28;
    const auto t157 = t156 * t34;
    const auto t158 = t0 + t1 - t17 - t22 - t25 + t7 + t84;
    const auto t159 = t28 * (-t115 * t155 + t155 * t157 * t41 + t158 + t35);
    const auto t160 = t0_y - t2_y;
    const auto t161 = t160 * t29;
    const auto t162 = t17 * t1_y;
    const auto t163 = t1_y * t25;
    const auto t164 = t17 * t2_y;
    const auto t165 = t25 * t2_y;
    const auto t166 = t0_y * t4;
    const auto t167 = t0_y * t8;
    const auto t168 = t0_y * t1;
    const auto t169 = t0_y * t7;
    const auto t170 = t0_y * t12 + t0_y * t15;
    const auto t171 = -t103 - t105 + t107 + t109 + t110 - 2 * t111 - 2 * t112
        + t162 + t163 - t164 - t165 + t166 + t167 - t168 - t169 + t170 - t96
        - t97 + t98 + t99;
    const auto t172 =
        t28 * (-t115 * t171 - t161 + 2 * t171 * t28 * t34 * t41 + t93);
    const auto t173 = t0_z - t2_z;
    const auto t174 = t173 * t29;
    const auto t175 = t17 * t1_z;
    const auto t176 = t1_z * t22;
    const auto t177 = t17 * t2_z;
    const auto t178 = t22 * t2_z;
    const auto t179 = t0_z * t4;
    const auto t180 = t0_z * t3;
    const auto t181 = t0_z * t1;
    const auto t182 = t0 * t0_z;
    const auto t183 = t0_z * t11 + t0_z * t15;
    const auto t184 = -t126 - t127 + t128 + t129 - t133 - t135 + t137 + t139
        + t140 - 2 * t141 - 2 * t142 + t175 + t176 - t177 - t178 + t179 + t180
        - t181 - t182 + t183;
    const auto t185 =
        t28 * (-t115 * t184 + t123 - t174 + 2 * t184 * t28 * t34 * t41);
    const auto t186 = -t146 - t147 + t148 + t149 - t150 - t151 + t152 + t153
        + t154 - t42 - t43 + t48 + t49 - t50 - t51 + t54 + t56 - 2 * t58
        - 2 * t59 + t63;
    const auto t187 = t186 * t41;
    const auto t188 =
        t28 * (-t115 * t186 + t157 * t187 - std::pow(t29, 2) + t41);
    const auto t189 = t29 * t91;
    const auto t190 = t100 + t101 - t102 - t104 + t106 + t108 + t113 - t162
        - t163 + t164 + t165 - t166 - t167 + t168 + t169 + t170 - t2 * t2_y
        - t2_y * t9 - t94 - t95;
    const auto t191 = t28 * (-t115 * t190 - t189 + 2 * t190 * t28 * t34 * t41);
    const auto t192 = t121 * t29;
    const auto t193 = -t124 - t125 + t130 + t131 - t132 - t134 + t136 + t138
        + t143 - t175 - t176 + t177 + t178 - t179 - t180 + t181 + t182 + t183
        - t2 * t2_z - t2_z * t5;
    const auto t194 = t28 * (-t115 * t193 - t192 + 2 * t193 * t28 * t34 * t41);
    const auto t195 = 2 * t161;
    const auto t196 = t156 * t91;
    const auto t197 =
        t28 * (2 * t160 * t28 * t41 * t64 - t195 - t196 * t77 - t33 * t91);
    const auto t198 = t160 * t91;
    const auto t199 = t196 * t76;
    const auto t200 = t28
        * (2 * t114 * t160 * t28 * t41 - t114 * t199 - 2 * t198 - t85
           - t90 * t91);
    const auto t201 = t121 * t160;
    const auto t202 = 2 * t201;
    const auto t203 =
        t28 * (-t120 * t91 + 2 * t144 * t160 * t28 * t41 - t144 * t199 - t202);
    const auto t204 = t156 * t160;
    const auto t205 = t204 * t41;
    const auto t206 = t28 * (-t155 * t199 + t155 * t205 + t195 - t92);
    const auto t207 = t28 * (t158 - t171 * t199 + t171 * t205 + t198);
    const auto t208 = t173 * t91;
    const auto t209 =
        t28 * (2 * t160 * t184 * t28 * t41 - t184 * t199 + t202 - t208);
    const auto t210 = t28 * (2 * t160 * t186 * t28 * t41 - t186 * t199 - t189);
    const auto t211 =
        t28 * (-t190 * t199 + t190 * t205 + t41 - std::pow(t91, 2));
    const auto t212 = t121 * t91;
    const auto t213 = t28 * (2 * t160 * t193 * t28 * t41 - t193 * t199 - t212);
    const auto t214 = 2 * t174;
    const auto t215 = t121 * t156;
    const auto t216 =
        t28 * (-t121 * t33 + 2 * t173 * t28 * t41 * t64 - t214 - t215 * t77);
    const auto t217 = 2 * t208;
    const auto t218 = t215 * t76;
    const auto t219 =
        t28 * (2 * t114 * t173 * t28 * t41 - t114 * t218 - t121 * t90 - t217);
    const auto t220 = t121 * t173;
    const auto t221 = t28
        * (-t120 * t121 + 2 * t144 * t173 * t28 * t41 - t144 * t218 - 2 * t220
           - t85);
    const auto t222 = t156 * t173;
    const auto t223 = t222 * t41;
    const auto t224 = t28 * (-t122 - t155 * t218 + t155 * t223 + t214);
    const auto t225 = t28 * (-t171 * t218 + t171 * t223 - t201 + t217);
    const auto t226 = t28 * (t158 - t184 * t218 + t184 * t223 + t220);
    const auto t227 = t28 * (2 * t173 * t186 * t28 * t41 - t186 * t218 - t192);
    const auto t228 = t28 * (2 * t173 * t190 * t28 * t41 - t190 * t218 - t212);
    const auto t229 =
        t28 * (-std::pow(t121, 2) - t193 * t218 + t193 * t223 + t41);
    const auto t230 = p_x + t32;
    const auto t231 = p_x - t0_x;
    const auto t232 = t231 * t29;
    const auto t233 = p_y - t0_y;
    const auto t234 = t233 * t91;
    const auto t235 = p_z - t0_z;
    const auto t236 = t121 * t235;
    const auto t237 = t234 + t236;
    const auto t238 = t232 + t237;
    const auto t239 = t238 * t76;
    const auto t240 = -t13;
    const auto t241 = -t14;
    const auto t242 = t20 + t21;
    const auto t243 = t240 + t241 + t242 + t36;
    const auto t244 = t231 * t34;
    const auto t245 = t160 * t233;
    const auto t246 = t173 * t235;
    const auto t247 = t245 + t246;
    const auto t248 = t244 + t247;
    const auto t249 = std::pow(t64, 2);
    const auto t250 = std::pow(t27, -2);
    const auto t251 = p_x + t2_x + t31;
    const auto t252 = t156 * t230;
    const auto t253 = t248 * t41;
    const auto t254 = 4 * t250;
    const auto t255 = 4 * t28;
    const auto t256 = t248 * t255;
    const auto t257 = t256 * t29;
    const auto t258 = -t232;
    const auto t259 = -t234;
    const auto t260 = -t236;
    const auto t261 = t258 + t259 + t260;
    const auto t262 = t251 * t78 - t257 * t64 + t261;
    const auto t263 = t68 + t7;
    const auto t264 = t0 + t70;
    const auto t265 = t263 + t264;
    const auto t266 = t12 - t19;
    const auto t267 = t11 - t18;
    const auto t268 = t265 + t266 + t267;
    const auto t269 = t1 + t73;
    const auto t270 = t15 - t23;
    const auto t271 = t248 + t268 + t269 + t270;
    const auto t272 = p_y + t89;
    const auto t273 = p_y + t2_y + t88;
    const auto t274 = t1_x * t1_y;
    const auto t275 = t2_x * t2_y;
    const auto t276 = t1_x * t2_y;
    const auto t277 = -t276;
    const auto t278 = t1_y * t2_x;
    const auto t279 = -t278;
    const auto t280 = t274 + t275 + t277 + t279;
    const auto t281 = t252 * t76;
    const auto t282 = t156 * t77;
    const auto t283 = 8 * t250;
    const auto t284 = t283 * t64;
    const auto t285 = 2 * t251;
    const auto t286 = t256 * t64;
    const auto t287 = t285 * t91 - t286 * t91;
    const auto t288 = -t114 * t257 + t273 * t78;
    const auto t289 = t28
        * (8 * t114 * t238 * t250 * t64 * t76 + 2 * t114 * t238 * t28 * t33
           + 2 * t114 * t251 * t28 * t41 - t114 * t253 * t284 - t114 * t281
           - t156 * t253 * t280 - t230 * t90 + 2 * t238 * t28 * t280 * t76
           + 2 * t238 * t28 * t64 * t90 - t272 * t282 - t272 * t33
           + 2 * t273 * t28 * t41 * t64 - t287 - t288);
    const auto t290 = p_z + t119;
    const auto t291 = p_z + t118 + t2_z;
    const auto t292 = t1_x * t1_z;
    const auto t293 = t2_x * t2_z;
    const auto t294 = t1_x * t2_z;
    const auto t295 = -t294;
    const auto t296 = t1_z * t2_x;
    const auto t297 = -t296;
    const auto t298 = t292 + t293 + t295 + t297;
    const auto t299 = t121 * t285 - t121 * t286;
    const auto t300 = -t144 * t257 + t291 * t78;
    const auto t301 = t28
        * (-t120 * t230 + 2 * t120 * t238 * t28 * t64
           + 8 * t144 * t238 * t250 * t64 * t76 + 2 * t144 * t238 * t28 * t33
           + 2 * t144 * t251 * t28 * t41 - t144 * t253 * t284 - t144 * t281
           - t156 * t253 * t298 + 2 * t238 * t28 * t298 * t76
           + 2 * t28 * t291 * t41 * t64 - t282 * t290 - t290 * t33 - t299
           - t300);
    const auto t302 = t155 * t257;
    const auto t303 = t156 * (t242 + t265 + t83);
    const auto t304 = t238 * t64;
    const auto t305 = t251 * t41;
    const auto t306 = t155 * t156;
    const auto t307 = t238 * t33;
    const auto t308 = t155 * t284;
    const auto t309 = 2 * t244 + 2 * t245 + 2 * t246 + t76;
    const auto t310 = t28
        * (-t155 * t281 + t157 * t304 - t230 * t34 + t231 * t282 + t231 * t33
           + t239 * t303 + t239 * t308 - t253 * t303 - t253 * t308 + t262 + t302
           + t305 * t306 + t306 * t307 + t309);
    const auto t311 = t0_y * t1_x;
    const auto t312 = -t311;
    const auto t313 = t1_y * t30;
    const auto t314 = 2 * t278;
    const auto t315 = t276 - t314;
    const auto t316 = t0_y * t2_x;
    const auto t317 = t2_y * t30;
    const auto t318 = t316 - t317;
    const auto t319 = t156 * (t275 + t312 + t313 + t315 + t318);
    const auto t320 = t156 * t171;
    const auto t321 = t171 * t257;
    const auto t322 = t171 * t284;
    const auto t323 = t28
        * (-t160 * t230 - t171 * t281 + t204 * t304 + t233 * t282 + t233 * t33
           - t239 * t319 + t239 * t322 + t253 * t319 - t253 * t322 + t287
           + t305 * t320 + t307 * t320 + t321);
    const auto t324 = t0_z * t1_x;
    const auto t325 = -t324;
    const auto t326 = t1_z * t30;
    const auto t327 = 2 * t296;
    const auto t328 = t294 - t327;
    const auto t329 = t0_z * t2_x;
    const auto t330 = t2_z * t30;
    const auto t331 = t329 - t330;
    const auto t332 = t156 * (t293 + t325 + t326 + t328 + t331);
    const auto t333 = t156 * t184;
    const auto t334 = t184 * t257;
    const auto t335 = t184 * t284;
    const auto t336 = t28
        * (-t173 * t230 - t184 * t281 + t222 * t304 + t235 * t282 + t235 * t33
           - t239 * t332 + t239 * t335 + t253 * t332 - t253 * t335 + t299
           + t305 * t333 + t307 * t333 + t334);
    const auto t337 = t231 * t41;
    const auto t338 = t156 * t64;
    const auto t339 = -t186 * t257 + t237;
    const auto t340 = t28
        * (-t156 * t239 * t268 + 8 * t186 * t238 * t250 * t64 * t76
           + 2 * t186 * t238 * t28 * t33 + 2 * t186 * t251 * t28 * t41
           - t186 * t253 * t284 - t186 * t281 - t230 * t29
           + 2 * t238 * t28 * t29 * t64 + 2 * t248 * t268 * t28 * t41 - t258
           - t337 * t338 - t339 - t41);
    const auto t341 = t233 * t29;
    const auto t342 = 2 * t341;
    const auto t343 = -t316;
    const auto t344 = 2 * t276;
    const auto t345 = t278 - t344;
    const auto t346 = t311 - t313;
    const auto t347 = t156 * (t274 + t317 + t343 + t345 + t346);
    const auto t348 = t233 * t41;
    const auto t349 = t156 * t305;
    const auto t350 = t190 * t257;
    const auto t351 = t190 * t238;
    const auto t352 = t156 * t33;
    const auto t353 = t190 * t284;
    const auto t354 = t28
        * (-t190 * t281 + t190 * t349 + t196 * t304 - t230 * t91 - t239 * t347
           + t239 * t353 + t253 * t347 - t253 * t353 - t338 * t348 + t342 + t350
           + t351 * t352);
    const auto t355 = t235 * t29;
    const auto t356 = 2 * t355;
    const auto t357 = -t329;
    const auto t358 = 2 * t294;
    const auto t359 = t296 - t358;
    const auto t360 = t324 - t326;
    const auto t361 = t156 * (t292 + t330 + t357 + t359 + t360);
    const auto t362 = t235 * t41;
    const auto t363 = t193 * t257;
    const auto t364 = t193 * t238;
    const auto t365 = t193 * t284;
    const auto t366 = t28
        * (-t121 * t230 - t193 * t281 + t193 * t349 + t215 * t304 - t239 * t361
           + t239 * t365 + t253 * t361 - t253 * t365 - t338 * t362 + t352 * t364
           + t356 + t363);
    const auto t367 = -t16;
    const auto t368 = t21 + t24;
    const auto t369 = t19 + t23;
    const auto t370 = t241 + t367 + t368 + t369;
    const auto t371 = std::pow(t114, 2);
    const auto t372 = t114 * t156;
    const auto t373 = t272 * t76;
    const auto t374 = t261 + t271;
    const auto t375 = 2 * t91;
    const auto t376 = t256 * t91;
    const auto t377 = -t114 * t376 + t273 * t375;
    const auto t378 = t1_y * t1_z;
    const auto t379 = t2_y * t2_z;
    const auto t380 = t1_y * t2_z;
    const auto t381 = -t380;
    const auto t382 = t1_z * t2_y;
    const auto t383 = -t382;
    const auto t384 = t378 + t379 + t381 + t383;
    const auto t385 = t144 * t156;
    const auto t386 = t372 * t76;
    const auto t387 = t114 * t283;
    const auto t388 = 2 * t121;
    const auto t389 = t121 * t256;
    const auto t390 = -t114 * t389 + t273 * t388;
    const auto t391 = -t144 * t376 + t291 * t375;
    const auto t392 = t28
        * (2 * t114 * t120 * t238 * t28 + 8 * t114 * t144 * t238 * t250 * t76
           + 2 * t114 * t28 * t291 * t41 - t120 * t272
           + 2 * t144 * t238 * t28 * t90 - t144 * t253 * t387
           + 2 * t144 * t273 * t28 * t41 - t156 * t253 * t384
           + 2 * t238 * t28 * t384 * t76 - t290 * t386 - t290 * t90
           - t373 * t385 - t390 - t391);
    const auto t393 = t0_x * t1_y;
    const auto t394 = -t393;
    const auto t395 = t1_x * t87;
    const auto t396 = t0_x * t2_y;
    const auto t397 = t2_x * t87;
    const auto t398 = t396 - t397;
    const auto t399 = t156 * (t275 + t345 + t394 + t395 + t398);
    const auto t400 = t114 * t238;
    const auto t401 = t273 * t41;
    const auto t402 = t238 * t90;
    const auto t403 = t155 * t376;
    const auto t404 = t155 * t387;
    const auto t405 = t28
        * (t157 * t400 + t231 * t386 + t231 * t90 - t239 * t399 + t239 * t404
           + t253 * t399 - t253 * t404 - t272 * t34 + t288 - t306 * t373
           + t306 * t401 + t306 * t402 + t403);
    const auto t406 = t171 * t376;
    const auto t407 = t263 + t269;
    const auto t408 = t156 * (t368 + t407 + t80 + t82);
    const auto t409 = t171 * t387;
    const auto t410 = t261 + t309;
    const auto t411 = t28
        * (-t160 * t272 + t204 * t400 + t233 * t386 + t233 * t90 + t239 * t408
           + t239 * t409 - t253 * t408 - t253 * t409 - t320 * t373 + t320 * t401
           + t320 * t402 + t377 + t406 + t410);
    const auto t412 = t0_z * t1_y;
    const auto t413 = -t412;
    const auto t414 = t1_z * t87;
    const auto t415 = 2 * t382;
    const auto t416 = t380 - t415;
    const auto t417 = t0_z * t2_y;
    const auto t418 = t2_z * t87;
    const auto t419 = t417 - t418;
    const auto t420 = t156 * (t379 + t413 + t414 + t416 + t419);
    const auto t421 = t184 * t376;
    const auto t422 = t184 * t387;
    const auto t423 = t28
        * (-t173 * t272 + t222 * t400 + t235 * t386 + t235 * t90 - t239 * t420
           + t239 * t422 + t253 * t420 - t253 * t422 - t333 * t373 + t333 * t401
           + t333 * t402 + t390 + t421);
    const auto t424 = t231 * t91;
    const auto t425 = 2 * t424;
    const auto t426 = -t396;
    const auto t427 = t393 - t395;
    const auto t428 = t156 * (t274 + t315 + t397 + t426 + t427);
    const auto t429 = t156 * t187;
    const auto t430 = t186 * t376;
    const auto t431 = t186 * t238;
    const auto t432 = t156 * t90;
    const auto t433 = t156 * t373;
    const auto t434 = t186 * t387;
    const auto t435 = t28
        * (-t186 * t433 - t239 * t428 + t239 * t434 + t253 * t428 - t253 * t434
           - t272 * t29 + t273 * t429 - t337 * t372 + t400 * t79 + t425 + t430
           + t431 * t432);
    const auto t436 = t266 + t270 + t407;
    const auto t437 = -t190 * t376 + t232 + t236;
    const auto t438 = t28
        * (8 * t114 * t190 * t238 * t250 * t76 + 2 * t114 * t238 * t28 * t91
           - t156 * t239 * t436 + 2 * t190 * t238 * t28 * t90
           - t190 * t253 * t387 + 2 * t190 * t273 * t28 * t41 - t190 * t433
           + 2 * t248 * t28 * t41 * t436 - t259 - t272 * t91 - t348 * t372 - t41
           - t437);
    const auto t439 = t235 * t91;
    const auto t440 = 2 * t439;
    const auto t441 = -t417;
    const auto t442 = 2 * t380;
    const auto t443 = t382 - t442;
    const auto t444 = t412 - t414;
    const auto t445 = t156 * (t378 + t418 + t441 + t443 + t444);
    const auto t446 = t193 * t376;
    const auto t447 = t193 * t387;
    const auto t448 = t28
        * (-t121 * t272 + t156 * t193 * t401 - t193 * t433 + t215 * t400
           - t239 * t445 + t239 * t447 + t253 * t445 - t253 * t447 - t362 * t372
           + t364 * t432 + t440 + t446);
    const auto t449 = t20 + t24;
    const auto t450 = t18 + t23;
    const auto t451 = t240 + t367 + t449 + t450;
    const auto t452 = std::pow(t144, 2);
    const auto t453 = t290 * t76;
    const auto t454 = -t144 * t389 + t291 * t388;
    const auto t455 = t0_x * t1_z;
    const auto t456 = -t455;
    const auto t457 = t117 * t1_x;
    const auto t458 = t0_x * t2_z;
    const auto t459 = t117 * t2_x;
    const auto t460 = t458 - t459;
    const auto t461 = t156 * (t293 + t359 + t456 + t457 + t460);
    const auto t462 = t385 * t76;
    const auto t463 = t144 * t238;
    const auto t464 = t291 * t41;
    const auto t465 = t120 * t238;
    const auto t466 = t155 * t389;
    const auto t467 = t144 * t283;
    const auto t468 = t155 * t467;
    const auto t469 = t28
        * (t120 * t231 + t157 * t463 + t231 * t462 - t239 * t461 + t239 * t468
           + t253 * t461 - t253 * t468 - t290 * t34 + t300 - t306 * t453
           + t306 * t464 + t306 * t465 + t466);
    const auto t470 = t0_y * t1_z;
    const auto t471 = -t470;
    const auto t472 = t117 * t1_y;
    const auto t473 = t0_y * t2_z;
    const auto t474 = t117 * t2_y;
    const auto t475 = t473 - t474;
    const auto t476 = t156 * (t379 + t443 + t471 + t472 + t475);
    const auto t477 = t171 * t389;
    const auto t478 = t171 * t467;
    const auto t479 = t28
        * (t120 * t233 - t160 * t290 + t204 * t463 + t233 * t462 - t239 * t476
           + t239 * t478 + t253 * t476 - t253 * t478 - t320 * t453 + t320 * t464
           + t320 * t465 + t391 + t477);
    const auto t480 = t184 * t389;
    const auto t481 = t264 + t269;
    const auto t482 = t156 * (t449 + t481 + t80 + t81);
    const auto t483 = t184 * t467;
    const auto t484 = t28
        * (t120 * t235 - t173 * t290 + t222 * t463 + t235 * t462 + t239 * t482
           + t239 * t483 - t253 * t482 - t253 * t483 - t333 * t453 + t333 * t464
           + t333 * t465 + t410 + t454 + t480);
    const auto t485 = t121 * t231;
    const auto t486 = 2 * t485;
    const auto t487 = -t458;
    const auto t488 = t455 - t457;
    const auto t489 = t156 * (t292 + t328 + t459 + t487 + t488);
    const auto t490 = t186 * t389;
    const auto t491 = t120 * t156;
    const auto t492 = t156 * t453;
    const auto t493 = t186 * t467;
    const auto t494 = t28
        * (-t186 * t492 - t239 * t489 + t239 * t493 + t253 * t489 - t253 * t493
           - t29 * t290 + t291 * t429 - t337 * t385 + t431 * t491 + t463 * t79
           + t486 + t490);
    const auto t495 = t121 * t233;
    const auto t496 = 2 * t495;
    const auto t497 = -t473;
    const auto t498 = t470 - t472;
    const auto t499 = t156 * (t378 + t416 + t474 + t497 + t498);
    const auto t500 = t190 * t389;
    const auto t501 = t190 * t467;
    const auto t502 = t28
        * (t156 * t190 * t464 - t190 * t492 + t196 * t463 - t239 * t499
           + t239 * t501 + t253 * t499 - t253 * t501 - t290 * t91 - t348 * t385
           + t351 * t491 + t496 + t500);
    const auto t503 = t267 + t270 + t481;
    const auto t504 = -t193 * t389 + t232 + t234;
    const auto t505 = t28
        * (2 * t120 * t193 * t238 * t28 + 2 * t121 * t144 * t238 * t28
           - t121 * t290 + 8 * t144 * t193 * t238 * t250 * t76
           - t156 * t239 * t503 - t193 * t253 * t467
           + 2 * t193 * t28 * t291 * t41 - t193 * t492
           + 2 * t248 * t28 * t41 * t503 - t260 - t362 * t385 - t41 - t504);
    const auto t506 = -t10 + t25;
    const auto t507 = t22 - t6;
    const auto t508 = t242 + t506 + t507;
    const auto t509 = std::pow(t155, 2);
    const auto t510 = t0_x * t0_y;
    const auto t511 = t156 * (t275 + t343 + t426 + t510);
    const auto t512 = t157 * t238;
    const auto t513 = t155 * t238;
    const auto t514 = t231 * t76;
    const auto t515 = t306 * t76;
    const auto t516 = t155 * t283;
    const auto t517 = t171 * t516;
    const auto t518 = t28
        * (t160 * t231 + t171 * t512 + t204 * t513 + t233 * t34 + t233 * t515
           + t239 * t511 + t239 * t517 - t253 * t511 - t253 * t517 + t320 * t514
           - t321 - t403);
    const auto t519 = t0_x * t0_z;
    const auto t520 = t156 * (t293 + t357 + t487 + t519);
    const auto t521 = t184 * t516;
    const auto t522 = t28
        * (t173 * t231 + t184 * t512 + t222 * t513 + t235 * t34 + t235 * t515
           + t239 * t520 + t239 * t521 - t253 * t520 - t253 * t521 + t333 * t514
           - t334 - t466);
    const auto t523 = t156 * t72;
    const auto t524 = t156 * t514;
    const auto t525 = t186 * t516;
    const auto t526 = t28
        * (t157 * t431 + t186 * t524 + t239 * t523 + t239 * t525 - t253 * t523
           - t253 * t525 - t306 * t337 + t339 + t513 * t79);
    const auto t527 = t156 * (t279 + t318 + t344 + t427 + t510);
    const auto t528 = t190 * t516;
    const auto t529 = t28
        * (t157 * t351 + t190 * t524 + t196 * t513 - t239 * t527 + t239 * t528
           + t253 * t527 - t253 * t528 - t306 * t348 - t342 - t350 + t424);
    const auto t530 = t156 * (t297 + t331 + t358 + t488 + t519);
    const auto t531 = t193 * t516;
    const auto t532 = t28
        * (t157 * t364 + t193 * t524 + t215 * t513 - t239 * t530 + t239 * t531
           + t253 * t530 - t253 * t531 - t306 * t362 - t356 - t363 + t485);
    const auto t533 = t17 - t26;
    const auto t534 = t368 + t506 + t533;
    const auto t535 = std::pow(t171, 2);
    const auto t536 = t0_y * t0_z;
    const auto t537 = t156 * (t379 + t441 + t497 + t536);
    const auto t538 = t171 * t238;
    const auto t539 = t233 * t76;
    const auto t540 = t235 * t76;
    const auto t541 = t171 * t283;
    const auto t542 = t184 * t541;
    const auto t543 = t28
        * (t160 * t235 + t173 * t233 + t184 * t204 * t238 + t222 * t538
           + t239 * t537 + t239 * t542 - t253 * t537 - t253 * t542 + t320 * t540
           + t333 * t539 - t421 - t477);
    const auto t544 = t156 * (t277 + t314 + t346 + t398 + t510);
    const auto t545 = t156 * t539;
    const auto t546 = t186 * t541;
    const auto t547 = t28
        * (t186 * t545 + t204 * t431 - t239 * t544 + t239 * t546 + t253 * t544
           - t253 * t546 - t320 * t337 + t341 - t425 - t430 + t538 * t79);
    const auto t548 = t156 * (t66 + t69 + t75);
    const auto t549 = t190 * t541;
    const auto t550 = t28
        * (t190 * t545 + t196 * t538 + t204 * t351 + t239 * t548 + t239 * t549
           - t253 * t548 - t253 * t549 - t320 * t348 + t437);
    const auto t551 = t156 * (t383 + t419 + t442 + t498 + t536);
    const auto t552 = t193 * t541;
    const auto t553 = t28
        * (t193 * t545 + t204 * t364 + t215 * t538 - t239 * t551 + t239 * t552
           + t253 * t551 - t253 * t552 - t320 * t362 - t440 - t446 + t495);
    const auto t554 = t449 + t507 + t533;
    const auto t555 = std::pow(t184, 2);
    const auto t556 = t156 * (t295 + t327 + t360 + t460 + t519);
    const auto t557 = t184 * t238;
    const auto t558 = t156 * t540;
    const auto t559 = t184 * t283;
    const auto t560 = t186 * t253;
    const auto t561 = t239 * t559;
    const auto t562 = t28
        * (t186 * t558 + t186 * t561 + t222 * t431 - t239 * t556 + t253 * t556
           - t333 * t337 + t355 - t486 - t490 + t557 * t79 - t559 * t560);
    const auto t563 = t156 * (t381 + t415 + t444 + t475 + t536);
    const auto t564 = t253 * t559;
    const auto t565 = t28
        * (t190 * t558 + t190 * t561 - t190 * t564 + t196 * t557 + t222 * t351
           - t239 * t563 + t253 * t563 - t333 * t348 + t439 - t496 - t500);
    const auto t566 = t156 * (t65 + t71 + t75);
    const auto t567 = t28
        * (t193 * t558 + t193 * t561 - t193 * t564 + t215 * t557 + t222 * t364
           + t239 * t566 - t253 * t566 - t333 * t362 + t504);
    const auto t568 = std::pow(t186, 2);
    const auto t569 = 2 * t250;
    const auto t570 = t274 + t312 + t394 + t510;
    const auto t571 = t190 * t255;
    const auto t572 = t569
        * (4 * t186 * t190 * t238 * t28 * t76 + t186 * t238 * t91 - t187 * t233
           + t190 * t238 * t29 - t190 * t337 + t238 * t570 * t76 - t253 * t570
           - t560 * t571);
    const auto t573 = t292 + t325 + t456 + t519;
    const auto t574 = t569
        * (t121 * t186 * t238 + 4 * t186 * t193 * t238 * t28 * t76 - t187 * t235
           + t193 * t238 * t29 - t193 * t255 * t560 - t193 * t337
           + t238 * t573 * t76 - t253 * t573);
    const auto t575 = t369 + t37 + t40;
    const auto t576 = std::pow(t190, 2);
    const auto t577 = t378 + t413 + t471 + t536;
    const auto t578 = t569
        * (t121 * t190 * t238 + 4 * t190 * t193 * t238 * t28 * t76 - t190 * t362
           + t193 * t238 * t91 - t193 * t253 * t571 - t193 * t348
           + t238 * t577 * t76 - t253 * t577);
    const auto t579 = t38 + t40 + t450;
    const auto t580 = std::pow(t193, 2);
    hess[0] = 0;
    hess[1] = 0;
    hess[2] = 0;
    hess[3] = t86;
    hess[4] = t116;
    hess[5] = t145;
    hess[6] = t159;
    hess[7] = t172;
    hess[8] = t185;
    hess[9] = t188;
    hess[10] = t191;
    hess[11] = t194;
    hess[12] = 0;
    hess[13] = 0;
    hess[14] = 0;
    hess[15] = t197;
    hess[16] = t200;
    hess[17] = t203;
    hess[18] = t206;
    hess[19] = t207;
    hess[20] = t209;
    hess[21] = t210;
    hess[22] = t211;
    hess[23] = t213;
    hess[24] = 0;
    hess[25] = 0;
    hess[26] = 0;
    hess[27] = t216;
    hess[28] = t219;
    hess[29] = t221;
    hess[30] = t224;
    hess[31] = t225;
    hess[32] = t226;
    hess[33] = t227;
    hess[34] = t228;
    hess[35] = t229;
    hess[36] = t86;
    hess[37] = t197;
    hess[38] = t216;
    hess[39] = t156
        * (-t230 * t33 + 4 * t238 * t249 * t250 * t76
           + 2 * t238 * t28 * t33 * t64 - t239 * t243 * t28
           + t243 * t248 * t28 * t41 - t249 * t253 * t254
           + 2 * t251 * t28 * t41 * t64 - t252 * t77 - t262 - t271);
    hess[40] = t289;
    hess[41] = t301;
    hess[42] = t310;
    hess[43] = t323;
    hess[44] = t336;
    hess[45] = t340;
    hess[46] = t354;
    hess[47] = t366;
    hess[48] = t116;
    hess[49] = t200;
    hess[50] = t219;
    hess[51] = t289;
    hess[52] = t156
        * (2 * t114 * t238 * t28 * t90 + 2 * t114 * t273 * t28 * t41
           + 4 * t238 * t250 * t371 * t76 - t239 * t28 * t370
           + t248 * t28 * t370 * t41 - t253 * t254 * t371 - t272 * t90
           - t372 * t373 - t374 - t377);
    hess[53] = t392;
    hess[54] = t405;
    hess[55] = t411;
    hess[56] = t423;
    hess[57] = t435;
    hess[58] = t438;
    hess[59] = t448;
    hess[60] = t145;
    hess[61] = t203;
    hess[62] = t221;
    hess[63] = t301;
    hess[64] = t392;
    hess[65] = t156
        * (2 * t120 * t144 * t238 * t28 - t120 * t290
           + 2 * t144 * t28 * t291 * t41 + 4 * t238 * t250 * t452 * t76
           - t239 * t28 * t451 + t248 * t28 * t41 * t451 - t253 * t254 * t452
           - t374 - t385 * t453 - t454);
    hess[66] = t469;
    hess[67] = t479;
    hess[68] = t484;
    hess[69] = t494;
    hess[70] = t502;
    hess[71] = t505;
    hess[72] = t159;
    hess[73] = t206;
    hess[74] = t224;
    hess[75] = t310;
    hess[76] = t405;
    hess[77] = t469;
    hess[78] = t156
        * (2 * t155 * t231 * t28 * t76 + 2 * t155 * t238 * t28 * t34
           + 4 * t238 * t250 * t509 * t76 - t239 * t28 * t508 - t247
           + t248 * t28 * t41 * t508 - t253 * t254 * t509 - t302);
    hess[79] = t518;
    hess[80] = t522;
    hess[81] = t526;
    hess[82] = t529;
    hess[83] = t532;
    hess[84] = t172;
    hess[85] = t207;
    hess[86] = t225;
    hess[87] = t323;
    hess[88] = t411;
    hess[89] = t479;
    hess[90] = t518;
    hess[91] = t156
        * (2 * t160 * t171 * t238 * t28 + 2 * t171 * t233 * t28 * t76
           + 4 * t238 * t250 * t535 * t76 - t239 * t28 * t534 - t244 - t246
           + t248 * t28 * t41 * t534 - t253 * t254 * t535 - t406);
    hess[92] = t543;
    hess[93] = t547;
    hess[94] = t550;
    hess[95] = t553;
    hess[96] = t185;
    hess[97] = t209;
    hess[98] = t226;
    hess[99] = t336;
    hess[100] = t423;
    hess[101] = t484;
    hess[102] = t522;
    hess[103] = t543;
    hess[104] = t156
        * (2 * t173 * t184 * t238 * t28 + 2 * t184 * t235 * t28 * t76
           + 4 * t238 * t250 * t555 * t76 - t239 * t28 * t554 - t244 - t245
           + t248 * t28 * t41 * t554 - t253 * t254 * t555 - t480);
    hess[105] = t562;
    hess[106] = t565;
    hess[107] = t567;
    hess[108] = t188;
    hess[109] = t210;
    hess[110] = t227;
    hess[111] = t340;
    hess[112] = t435;
    hess[113] = t494;
    hess[114] = t526;
    hess[115] = t547;
    hess[116] = t562;
    hess[117] = t569
        * (2 * t186 * t238 * t29 - 2 * t186 * t337 + 4 * t238 * t28 * t568 * t76
           - t239 * t39 + t248 * t39 * t41 - t253 * t255 * t568);
    hess[118] = t572;
    hess[119] = t574;
    hess[120] = t191;
    hess[121] = t211;
    hess[122] = t228;
    hess[123] = t354;
    hess[124] = t438;
    hess[125] = t502;
    hess[126] = t529;
    hess[127] = t550;
    hess[128] = t565;
    hess[129] = t572;
    hess[130] = t569
        * (2 * t190 * t238 * t91 - 2 * t190 * t348 + 4 * t238 * t28 * t576 * t76
           - t239 * t575 + t248 * t41 * t575 - t253 * t255 * t576);
    hess[131] = t578;
    hess[132] = t194;
    hess[133] = t213;
    hess[134] = t229;
    hess[135] = t366;
    hess[136] = t448;
    hess[137] = t505;
    hess[138] = t532;
    hess[139] = t553;
    hess[140] = t567;
    hess[141] = t574;
    hess[142] = t578;
    hess[143] = t569
        * (2 * t121 * t193 * t238 - 2 * t193 * t362
           + 4 * t238 * t28 * t580 * t76 - t239 * t579 + t248 * t41 * t579
           - t253 * t255 * t580);
}

void face_term_aux_fast_gradient(
    double t0_x,
    double t0_y,
    double t0_z,
    double t1_x,
    double t1_y,
    double t1_z,
    double t2_x,
    double t2_y,
    double t2_z,
    double p_x,
    double p_y,
    double p_z,
    double d,
    double grad[13])
{
    const auto t0 = t0_x - t1_x;
    const auto t1 = -t2_z;
    const auto t2 = t0_z + t1;
    const auto t3 = -t2_x;
    const auto t4 = t0_x + t3;
    const auto t5 = t0_z - t1_z;
    const auto t6 = t0 * t2 - t4 * t5;
    const auto t7 = -t2_y;
    const auto t8 = t0_y + t7;
    const auto t9 = t0 * t8;
    const auto t10 = t0_y - t1_y;
    const auto t11 = t10 * t4;
    const auto t12 = -t11 + t9;
    const auto t13 = t10 * t2 - t5 * t8;
    const auto t14 = std::pow(t12, 2) + std::pow(t13, 2);
    const auto t15 = t14 + std::pow(t6, 2);
    const auto t16 = std::pow(d * t15, -1.0 / 2.0);
    const auto t17 = p_y - t0_y;
    const auto t18 = t1 + t1_z;
    const auto t19 = p_z - t0_z;
    const auto t20 = t1_y + t7;
    const auto t21 = p_x - t0_x;
    const auto t22 = t12 * t19 + t13 * t21;
    const auto t23 = (-t17 * t6 + t22) / t15;
    const auto t24 = t1_x + t3;
    const auto t25 = t4 * t5;
    const auto t26 = t0 * t2 - t25;
    const auto t27 = t14 + std::pow(t26, 2);
    const auto t28 = std::pow(d * t27, -1.0 / 2.0);
    const auto t29 = (-t17 * t26 + t22) / t27;
    grad[0] =
        t16 * (-t13 - t17 * t18 + t19 * t20 - t23 * (t12 * t20 + t18 * t6));
    grad[1] = t16
        * (t0 * t2 + t18 * t21 - t19 * t24 - t23 * (-t12 * t24 + t13 * t18)
           - t25);
    grad[2] =
        t16 * (t11 + t17 * t24 - t20 * t21 + t23 * (t13 * t20 + t24 * t6) - t9);
    grad[3] = t28 * (t17 * t2 - t19 * t8 + t29 * (t12 * t8 + t2 * t6));
    grad[4] = -t28 * (-t19 * t4 + t2 * t21 + t29 * (t12 * t4 - t13 * t2));
    grad[5] = t28 * (-t17 * t4 + t21 * t8 - t29 * (t13 * t8 + t26 * t4));
    grad[6] = -t28 * (-t10 * t19 + t17 * t5 + t29 * (t10 * t12 + t26 * t5));
    grad[7] = t28 * (-t0 * t19 + t21 * t5 + t29 * (t0 * t12 - t13 * t5));
    grad[8] = t28 * (t0 * t17 - t10 * t21 + t29 * (t0 * t6 + t10 * t13));
    grad[9] = t13 * t28;
    grad[10] = -t26 * t28;
    grad[11] = t12 * t28;
    grad[12] = -1.0 / 2.0 * t15 * t28 * t29 / d;
}

// hess is (169×1) flattened in column-major order
void face_term_aux_fast_hessian(
    double t0_x,
    double t0_y,
    double t0_z,
    double t1_x,
    double t1_y,
    double t1_z,
    double t2_x,
    double t2_y,
    double t2_z,
    double p_x,
    double p_y,
    double p_z,
    double d,
    double hess[169])
{
    const auto t0 = -t2_y;
    const auto t1 = t0 + t1_y;
    const auto t2 = std::pow(t1, 2);
    const auto t3 = -t2_z;
    const auto t4 = t1_z + t3;
    const auto t5 = std::pow(t4, 2);
    const auto t6 = -t0_y;
    const auto t7 = p_y + t6;
    const auto t8 = -t1_x;
    const auto t9 = t0_x + t8;
    const auto t10 = -t9;
    const auto t11 = t0_z + t3;
    const auto t12 = -t11;
    const auto t13 = -t2_x;
    const auto t14 = t0_x + t13;
    const auto t15 = -t14;
    const auto t16 = -t1_z;
    const auto t17 = t0_z + t16;
    const auto t18 = -t17;
    const auto t19 = t10 * t12 - t15 * t18;
    const auto t20 = -t0_x;
    const auto t21 = p_x + t20;
    const auto t22 = -t1_y;
    const auto t23 = t0_y + t22;
    const auto t24 = t11 * t23;
    const auto t25 = t0 + t0_y;
    const auto t26 = t17 * t25;
    const auto t27 = -t26;
    const auto t28 = t24 + t27;
    const auto t29 = -t0_z;
    const auto t30 = p_z + t29;
    const auto t31 = t25 * t9;
    const auto t32 = t14 * t23;
    const auto t33 = -t32;
    const auto t34 = t31 + t33;
    const auto t35 = t21 * t28 + t30 * t34;
    const auto t36 = -t19 * t7 + t35;
    const auto t37 = -t1 * t30 + t28 + t4 * t7;
    const auto t38 = -t37;
    const auto t39 = t1 * t34;
    const auto t40 = t19 * t4 + t39;
    const auto t41 = 2 * t40;
    const auto t42 = std::pow(t28, 2) + std::pow(t34, 2);
    const auto t43 = std::pow(t19, 2) + t42;
    const auto t44 = 1.0 / t43;
    const auto t45 = t36 * t44;
    const auto t46 = -t25;
    const auto t47 = -t23;
    const auto t48 = t10 * t46 - t15 * t47;
    const auto t49 = -t19;
    const auto t50 = t1 * t48 - t4 * t49;
    const auto t51 = t36 * t50;
    const auto t52 = t44 * t51;
    const auto t53 = std::pow(d * t43, -1.0 / 2.0);
    const auto t54 = t44 * t53;
    const auto t55 = t13 + t1_x;
    const auto t56 = -t28 * t4 + t34 * t55;
    const auto t57 = -t56;
    const auto t58 = 2 * t57;
    const auto t59 = t21 * t4;
    const auto t60 = t30 * t55;
    const auto t61 = t11 * t9;
    const auto t62 = t14 * t17;
    const auto t63 = -t62;
    const auto t64 = t61 + t63;
    const auto t65 = t59 - t60 + t64;
    const auto t66 = t36 * t55;
    const auto t67 = t40 * t45;
    const auto t68 = t1 * t66 - t38 * t57 - t40 * t65 + t57 * t67;
    const auto t69 = t1 * t21;
    const auto t70 = -t34 + t55 * t7 - t69;
    const auto t71 = t40 * t70;
    const auto t72 = t1 * t28 + t19 * t55;
    const auto t73 = t67 * t72;
    const auto t74 = 2 * t72;
    const auto t75 = t11 * t7 - t25 * t30;
    const auto t76 = t40 * t75;
    const auto t77 = t11 * t19 + t25 * t34;
    const auto t78 = t1 * t25;
    const auto t79 = t11 * t4;
    const auto t80 = t78 + t79;
    const auto t81 = t67 * t77;
    const auto t82 = 2 * t77;
    const auto t83 = -t11 * t28 + t14 * t34;
    const auto t84 = std::pow(t43, -2);
    const auto t85 = -p_z;
    const auto t86 = t11 * t21 - t14 * t30;
    const auto t87 = -t86;
    const auto t88 = t40 * t44;
    const auto t89 = t38 * t44;
    const auto t90 = t2_z - t36 * t40 * t83 * t84 + t45 * (t1 * t14 + t34)
        + t83 * t89 + t85 + t87 * t88;
    const auto t91 = t25 * t28;
    const auto t92 = t14 * t19 + t91;
    const auto t93 = t21 * t25;
    const auto t94 = t14 * t7;
    const auto t95 = t93 - t94;
    const auto t96 = p_y + t0;
    const auto t97 = -t36 * t40 * t84 * t92 + t45 * (t14 * t4 + t64) + t88 * t95
        + t89 * t92 + t96;
    const auto t98 = t23 * t34;
    const auto t99 = t17 * t19 + t98;
    const auto t100 = t17 * t7 - t23 * t30;
    const auto t101 = -t100;
    const auto t102 = t1 * t23;
    const auto t103 = t17 * t4;
    const auto t104 =
        t101 * t40 - t36 * t40 * t44 * t99 + t36 * (t102 + t103) + t38 * t99;
    const auto t105 = -t17 * t28 + t34 * t9;
    const auto t106 = -t105;
    const auto t107 = t17 * t21;
    const auto t108 = t30 * t9;
    const auto t109 = t107 - t108;
    const auto t110 = p_z - t106 * t36 * t40 * t84 + t106 * t89 + t109 * t88
        + t16 + t45 * (-t1 * t9 - t34);
    const auto t111 = t19 * t9 + t23 * t28;
    const auto t112 = 2 * t111;
    const auto t113 = -p_y;
    const auto t114 = t21 * t23;
    const auto t115 = -t114 + t7 * t9;
    const auto t116 = t36 * t84;
    const auto t117 = t111 * t116 * t40 - t111 * t38 * t44 + t113 + t115 * t88
        + t1_y + t45 * (-t4 * t9 - t64);
    const auto t118 = t39 + t4 * t64;
    const auto t119 = t42 + std::pow(t64, 2);
    const auto t120 = std::pow(d * t119, -1.0 / 2.0);
    const auto t121 = 1.0 / t119;
    const auto t122 = t121 * t28;
    const auto t123 = t120 * t122;
    const auto t124 = -t118 * t123;
    const auto t125 = t53 * (t16 + t19 * t88 + t2_z);
    const auto t126 = t121 * t34;
    const auto t127 = t22 + t2_y;
    const auto t128 = t120 * (-t118 * t126 - t127);
    const auto t129 = t12 * t47 - t18 * t46;
    const auto t130 = std::pow(t129, 2) + std::pow(t48, 2) + std::pow(t49, 2);
    const auto t131 = (1.0 / 2.0) * t130;
    const auto t132 = t130 * t44;
    const auto t133 = 1.0 / d;
    const auto t134 = t133 * t54;
    const auto t135 = t129 * t4 - t48 * t55;
    const auto t136 = t135 * t36;
    const auto t137 = t136 * t44;
    const auto t138 = std::pow(t55, 2);
    const auto t139 = t57 * t70;
    const auto t140 = t45 * t57;
    const auto t141 = t140 * t72;
    const auto t142 = t136 * t84;
    const auto t143 = t44 * t57;
    const auto t144 = -t31 + t32;
    const auto t145 = t116 * t57;
    const auto t146 = p_z + t143 * t75 + t145 * t77 + t3 - t44 * t65 * t77
        + t45 * (t144 + t25 * t55);
    const auto t147 = 2 * t83;
    const auto t148 = t14 * t55;
    const auto t149 = t140 * t83 + t36 * (t148 + t79) - t57 * t87 - t65 * t83;
    const auto t150 = -p_x;
    const auto t151 = t44 * t65;
    const auto t152 = t143 * t95 + t150 + t151 * t92 + t2_x
        - t36 * t57 * t84 * t92 + t45 * (t25 * t4 + t28);
    const auto t153 = t101 * t143 + t151 * t99 + t1_z - t36 * t57 * t84 * t99
        + t45 * (-t144 - t23 * t55) + t85;
    const auto t154 = t55 * t9;
    const auto t155 =
        -t106 * t36 * t44 * t57 + t106 * t65 + t109 * t57 + t36 * (t103 + t154);
    const auto t156 = p_x + t111 * t145 - t111 * t44 * t65 + t115 * t143
        + t45 * (-t23 * t4 - t28) + t8;
    const auto t157 = t120 * (t122 * t56 + t4);
    const auto t158 = t121 * t64;
    const auto t159 = t120 * t158;
    const auto t160 = -t159 * t56;
    const auto t161 = -t53 * (t143 * t34 + t55);
    const auto t162 = -t1 * t129 + t49 * t55;
    const auto t163 = t162 * t36;
    const auto t164 = t163 * t44;
    const auto t165 = t44 * t72;
    const auto t166 = t44 * t70;
    const auto t167 = -t61 + t62;
    const auto t168 = t116 * t72;
    const auto t169 =
        t165 * t75 + t166 * t77 + t168 * t77 - t45 * (t11 * t55 + t167) + t96;
    const auto t170 = -t24 + t26;
    const auto t171 = p_x + t13 + t166 * t83 + t168 * t83 - t44 * t72 * t87
        + t45 * (t1 * t11 + t170);
    const auto t172 = t70 * t92;
    const auto t173 = t148 + t78;
    const auto t174 = 2 * t92;
    const auto t175 = t45 * t72;
    const auto t176 = t175 * t92;
    const auto t177 = p_y - t101 * t44 * t72 + t166 * t99 + t168 * t99 + t22
        + t45 * (-t167 - t17 * t55);
    const auto t178 = t106 * t166 + t106 * t168 - t109 * t44 * t72 + t150 + t1_x
        + t45 * (-t1 * t17 - t170);
    const auto t179 =
        t111 * t175 + t111 * t70 + t115 * t72 - t36 * (t102 + t154);
    const auto t180 = t53 * (t127 + t165 * t28);
    const auto t181 = t120 * (-t158 * t72 - t2_x - t8);
    const auto t182 = t120 * t126;
    const auto t183 = t182 * t72;
    const auto t184 = t11 * t49 + t46 * t48;
    const auto t185 = t184 * t36;
    const auto t186 = t185 * t44;
    const auto t187 = t185 * t84;
    const auto t188 = std::pow(t25, 2);
    const auto t189 = std::pow(t11, 2);
    const auto t190 = t75 * t83;
    const auto t191 = t14 * t36;
    const auto t192 = t45 * t77;
    const auto t193 = t192 * t83;
    const auto t194 = t75 * t92;
    const auto t195 = t192 * t92;
    const auto t196 = t75 * t99;
    const auto t197 = t23 * t25;
    const auto t198 = t11 * t17;
    const auto t199 = t197 + t198;
    const auto t200 = 2 * t99;
    const auto t201 = t192 * t99;
    const auto t202 = t44 * t75;
    const auto t203 = t116 * t77;
    const auto t204 = t0_z + t106 * t202 + t106 * t203 - t109 * t44 * t77
        + t45 * (2 * t31 + t33) + t85;
    const auto t205 = t44 * t77;
    const auto t206 = t0_y + t113;
    const auto t207 =
        t111 * t202 + t111 * t203 + t115 * t205 + t206 - t45 * (2 * t61 + t63);
    const auto t208 = t123 * t77;
    const auto t209 = t120 * (-t158 * t77 - t29 - t2_z);
    const auto t210 = t2_y + t6;
    const auto t211 = t53 * (t205 * t34 + t210);
    const auto t212 = t12 * t129 + t14 * t48;
    const auto t213 = t212 * t36;
    const auto t214 = t213 * t44;
    const auto t215 = t213 * t84;
    const auto t216 = std::pow(t14, 2);
    const auto t217 = t45 * t83;
    const auto t218 = t11 * t25 * t36 + t217 * t92 - t83 * t95 - t87 * t92;
    const auto t219 = t44 * t83;
    const auto t220 = t101 * t219 + t30 - t36 * t83 * t84 * t99
        + t44 * t87 * t99 + t45 * (2 * t14 * t23 - t31);
    const auto t221 = 2 * t106;
    const auto t222 = t14 * t9;
    const auto t223 =
        t106 * t217 - t106 * t87 - t109 * t83 + t36 * (t198 + t222);
    const auto t224 = t0_x + t111 * t116 * t83 - t111 * t44 * t87 + t115 * t219
        + t150 + t45 * (2 * t24 + t27);
    const auto t225 = -t53 * (t11 + t219 * t28);
    const auto t226 = t159 * t83;
    const auto t227 = t20 + t2_x;
    const auto t228 = t120 * (-t126 * t83 - t227);
    const auto t229 = t129 * t25 + t15 * t49;
    const auto t230 = t229 * t36;
    const auto t231 = t230 * t44;
    const auto t232 = t44 * t95;
    const auto t233 = t44 * t92;
    const auto t234 = t101 * t233 + t206 + t232 * t99 - t36 * t84 * t92 * t99
        + t45 * (2 * t14 * t17 - t61);
    const auto t235 = t106 * t232 - t106 * t36 * t84 * t92 + t109 * t233 + t21
        + t45 * (2 * t17 * t25 - t24);
    const auto t236 = t115 * t92;
    const auto t237 = t197 + t222;
    const auto t238 = t45 * t92;
    const auto t239 = t111 * t238;
    const auto t240 = t14 * t64 + t91;
    const auto t241 = t120 * (-t122 * t240 - t210);
    const auto t242 = t53 * (t19 * t233 + t227);
    const auto t243 = -t182 * t240;
    const auto t244 = t18 * t49 + t23 * t48;
    const auto t245 = t244 * t36;
    const auto t246 = t245 * t44;
    const auto t247 = std::pow(t23, 2);
    const auto t248 = std::pow(t17, 2);
    const auto t249 = t36 * t9;
    const auto t250 = t45 * t99;
    const auto t251 = -t101 * t106 + t106 * t250 - t109 * t99 + t23 * t249;
    const auto t252 = t115 * t99;
    const auto t253 = t111 * t250;
    const auto t254 = t17 * t64 + t98;
    const auto t255 = -t123 * t254;
    const auto t256 = t53 * (t19 * t44 * t99 + t1_z + t29);
    const auto t257 = t1_y + t6;
    const auto t258 = t120 * (-t126 * t254 - t257);
    const auto t259 = t10 * t48 + t129 * t17;
    const auto t260 = t259 * t36;
    const auto t261 = t260 * t84;
    const auto t262 = t260 * t44;
    const auto t263 = std::pow(t9, 2);
    const auto t264 = t106 * t115;
    const auto t265 = t106 * t45;
    const auto t266 = t111 * t265;
    const auto t267 = t120 * (t105 * t122 + t17);
    const auto t268 = -t105 * t159;
    const auto t269 = -t53 * (t106 * t34 * t44 + t9);
    const auto t270 = t129 * t47 + t49 * t9;
    const auto t271 = t270 * t36;
    const auto t272 = t271 * t44;
    const auto t273 = t53 * (t111 * t28 * t44 + t257);
    const auto t274 = t120 * (-t111 * t158 - t1_x - t20);
    const auto t275 = t111 * t182;
    const auto t276 = t111 * t45;
    const auto t277 = t121 * t43;
    const auto t278 = t120 * t277;
    const auto t279 = (1.0 / 2.0) * t133 * t278;
    const auto t280 = -t279 * t28;
    const auto t281 = t279 * t64;
    const auto t282 = -t279 * t34;
    const auto t283 = t131 * t134;
    hess[0] = t54
        * (-t36 * (t2 + t5) - t38 * t41 + std::pow(t40, 2) * t45 + t41 * t52);
    hess[1] = t54 * (t52 * t58 + t68);
    hess[2] = t54 * (t36 * t4 * t55 + t38 * t72 - t52 * t74 - t71 - t73);
    hess[3] = t54 * (t36 * t80 + t38 * t77 - t52 * t82 - t76 - t81);
    hess[4] = t53 * (2 * t36 * t50 * t83 * t84 - t90);
    hess[5] = t53 * (2 * t36 * t50 * t84 * t92 - t97);
    hess[6] = t54 * (-t104 + 2 * t36 * t44 * t50 * t99);
    hess[7] = t53 * (2 * t106 * t36 * t50 * t84 - t110);
    hess[8] = t53 * (-t112 * t51 * t84 - t117);
    hess[9] = t124;
    hess[10] = t125;
    hess[11] = t128;
    hess[12] = t134 * (-t131 * t38 + t131 * t67 + t132 * t51 - t51);
    hess[13] = t54 * (t137 * t41 + t68);
    hess[14] = t54
        * (t137 * t58 - t36 * (t138 + t5) + t45 * std::pow(t57, 2) - t58 * t65);
    hess[15] = t54 * (t1 * t36 * t4 - t137 * t74 - t139 - t141 + t65 * t72);
    hess[16] = t53 * (-t142 * t82 - t146);
    hess[17] = t54 * (t137 * t147 + t149);
    hess[18] = t53 * (2 * t135 * t36 * t84 * t92 - t152);
    hess[19] = t53 * (2 * t135 * t36 * t84 * t99 - t153);
    hess[20] = t54 * (2 * t106 * t135 * t36 * t44 - t155);
    hess[21] = t53 * (-t112 * t142 - t156);
    hess[22] = t157;
    hess[23] = t160;
    hess[24] = t161;
    hess[25] = t134 * (t131 * t140 - t131 * t65 + t132 * t136 - t136);
    hess[26] = t54 * (t164 * t41 + t38 * t72 + t4 * t66 - t71 - t73);
    hess[27] = t54 * (t1 * t36 * t4 - t139 - t141 + t164 * t58 + t65 * t72);
    hess[28] = t54
        * (-t164 * t74 + t36 * t44 * std::pow(t72, 2) - t36 * (t138 + t2)
           + 2 * t70 * t72);
    hess[29] = t53 * (-t163 * t82 * t84 + t169);
    hess[30] = t53 * (2 * t162 * t36 * t83 * t84 - t171);
    hess[31] = t54 * (t164 * t174 - t172 + t173 * t36 - t176 + t72 * t95);
    hess[32] = t53 * (2 * t162 * t36 * t84 * t99 - t177);
    hess[33] = t53 * (2 * t106 * t162 * t36 * t84 - t178);
    hess[34] = t54 * (-t112 * t164 + t179);
    hess[35] = t180;
    hess[36] = t181;
    hess[37] = t183;
    hess[38] =
        t134 * (t130 * t162 * t36 * t44 - t131 * t175 - t131 * t70 - t163);
    hess[39] = t54 * (t186 * t41 + t36 * t80 + t38 * t77 - t76 - t81);
    hess[40] = t53 * (-t146 + 2 * t184 * t36 * t57 * t84);
    hess[41] = t53 * (t169 - t187 * t74);
    hess[42] = t54
        * (-t186 * t82 - t36 * (t188 + t189) + t45 * std::pow(t77, 2)
           + t75 * t82);
    hess[43] = t54 * (t147 * t186 - t190 + t191 * t25 - t193 + t77 * t87);
    hess[44] = t54 * (t11 * t191 + t174 * t186 - t194 - t195 + t77 * t95);
    hess[45] = t54 * (t101 * t77 + t186 * t200 - t196 + t199 * t36 - t201);
    hess[46] = t53 * (2 * t106 * t184 * t36 * t84 - t204);
    hess[47] = t53 * (-t112 * t187 + t207);
    hess[48] = t208;
    hess[49] = t209;
    hess[50] = t211;
    hess[51] =
        t134 * (t130 * t184 * t36 * t44 - t131 * t192 - t131 * t75 - t185);
    hess[52] = t53 * (2 * t212 * t36 * t40 * t84 - t90);
    hess[53] = t54 * (t149 + t214 * t58);
    hess[54] = t53 * (-t171 - t215 * t74);
    hess[55] = t54 * (t14 * t25 * t36 - t190 - t193 - t214 * t82 + t77 * t87);
    hess[56] = t54
        * (-t147 * t87 + 2 * t212 * t36 * t44 * t83
           + t36 * t44 * std::pow(t83, 2) - t36 * (t189 + t216));
    hess[57] = t54 * (t174 * t214 + t218);
    hess[58] = t53 * (2 * t212 * t36 * t84 * t99 - t220);
    hess[59] = t54 * (t214 * t221 + t223);
    hess[60] = t53 * (-t112 * t215 - t224);
    hess[61] = t225;
    hess[62] = t226;
    hess[63] = t228;
    hess[64] = t134
        * (t130 * t212 * t36 * t44 + (1.0 / 2.0) * t130 * t36 * t44 * t83
           - t131 * t87 - t213);
    hess[65] = t53 * (2 * t229 * t36 * t40 * t84 - t97);
    hess[66] = t53 * (-t152 + 2 * t229 * t36 * t57 * t84);
    hess[67] = t54 * (-t172 + t173 * t36 - t176 - t231 * t74 + t72 * t95);
    hess[68] = t54 * (t11 * t14 * t36 - t194 - t195 - t231 * t82 + t77 * t95);
    hess[69] = t54 * (t147 * t231 + t218);
    hess[70] = t54
        * (-t174 * t95 + 2 * t229 * t36 * t44 * t92
           + t36 * t44 * std::pow(t92, 2) - t36 * (t188 + t216));
    hess[71] = t53 * (2 * t229 * t36 * t84 * t99 - t234);
    hess[72] = t53 * (2 * t106 * t229 * t36 * t84 - t235);
    hess[73] = t54 * (t111 * t95 - t112 * t231 - t236 + t237 * t36 - t239);
    hess[74] = t241;
    hess[75] = t242;
    hess[76] = t243;
    hess[77] = t134
        * (t130 * t229 * t36 * t44 + (1.0 / 2.0) * t130 * t36 * t44 * t92
           - t131 * t95 - t230);
    hess[78] = t54 * (-t104 + 2 * t244 * t36 * t40 * t44);
    hess[79] = t53 * (-t153 + 2 * t244 * t36 * t57 * t84);
    hess[80] = t53 * (-t177 - t245 * t74 * t84);
    hess[81] = t54 * (t101 * t77 - t196 + t199 * t36 - t201 - t246 * t82);
    hess[82] = t53 * (-t220 + 2 * t244 * t36 * t83 * t84);
    hess[83] = t53 * (-t234 + 2 * t244 * t36 * t84 * t92);
    hess[84] = t54
        * (-t101 * t200 + 2 * t244 * t36 * t44 * t99
           + t36 * t44 * std::pow(t99, 2) - t36 * (t247 + t248));
    hess[85] = t54 * (t221 * t246 + t251);
    hess[86] = t54 * (t101 * t111 - t112 * t246 + t17 * t36 * t9 - t252 - t253);
    hess[87] = t255;
    hess[88] = t256;
    hess[89] = t258;
    hess[90] = t134
        * (-t101 * t131 + t130 * t244 * t36 * t44
           + (1.0 / 2.0) * t130 * t36 * t44 * t99 - t245);
    hess[91] = t53 * (-t110 + 2 * t259 * t36 * t40 * t84);
    hess[92] = t54 * (-t155 + 2 * t259 * t36 * t44 * t57);
    hess[93] = t53 * (-t178 - t261 * t74);
    hess[94] = t53 * (-t204 - t261 * t82);
    hess[95] = t54 * (t147 * t262 + t223);
    hess[96] = t53 * (-t235 + 2 * t259 * t36 * t84 * t92);
    hess[97] = t54 * (t200 * t262 + t251);
    hess[98] = t54
        * (std::pow(t106, 2) * t36 * t44 + 2 * t106 * t259 * t36 * t44
           - t109 * t221 - t36 * (t248 + t263));
    hess[99] =
        t54 * (t109 * t111 - t112 * t262 + t17 * t23 * t36 - t264 - t266);
    hess[100] = t267;
    hess[101] = t268;
    hess[102] = t269;
    hess[103] = t134
        * ((1.0 / 2.0) * t106 * t130 * t36 * t44 - t109 * t131
           + t130 * t259 * t36 * t44 - t260);
    hess[104] = t53 * (-t117 + 2 * t270 * t36 * t40 * t84);
    hess[105] = t53 * (-t156 + 2 * t270 * t36 * t57 * t84);
    hess[106] = t54 * (t179 - t272 * t74);
    hess[107] = t53 * (t207 - t271 * t82 * t84);
    hess[108] = t53 * (-t224 + 2 * t270 * t36 * t83 * t84);
    hess[109] = t54 * (t111 * t95 + t174 * t272 - t236 + t237 * t36 - t239);
    hess[110] = t54 * (t101 * t111 + t17 * t249 + t200 * t272 - t252 - t253);
    hess[111] =
        t54 * (t109 * t111 + t17 * t23 * t36 + t221 * t272 - t264 - t266);
    hess[112] = t54
        * (std::pow(t111, 2) * t45 + t112 * t115 - t112 * t272
           - t36 * (t247 + t263));
    hess[113] = t273;
    hess[114] = t274;
    hess[115] = t275;
    hess[116] =
        t134 * (-t115 * t131 + t130 * t270 * t36 * t44 - t131 * t276 - t271);
    hess[117] = t124;
    hess[118] = t157;
    hess[119] = t180;
    hess[120] = t208;
    hess[121] = t225;
    hess[122] = t241;
    hess[123] = t255;
    hess[124] = t267;
    hess[125] = t273;
    hess[126] = 0;
    hess[127] = 0;
    hess[128] = 0;
    hess[129] = t280;
    hess[130] = t125;
    hess[131] = t160;
    hess[132] = t181;
    hess[133] = t209;
    hess[134] = t226;
    hess[135] = t242;
    hess[136] = t256;
    hess[137] = t268;
    hess[138] = t274;
    hess[139] = 0;
    hess[140] = 0;
    hess[141] = 0;
    hess[142] = t281;
    hess[143] = t128;
    hess[144] = t161;
    hess[145] = t183;
    hess[146] = t211;
    hess[147] = t228;
    hess[148] = t243;
    hess[149] = t258;
    hess[150] = t269;
    hess[151] = t275;
    hess[152] = 0;
    hess[153] = 0;
    hess[154] = 0;
    hess[155] = t282;
    hess[156] = t283 * (t37 + t67);
    hess[157] = t283 * (t140 + t167 - t59 + t60);
    hess[158] = -t283 * (t144 + t175 + t55 * t7 - t69);
    hess[159] = -t283 * (t192 + t75);
    hess[160] = t283 * (t217 + t86);
    hess[161] = t283 * (t238 - t93 + t94);
    hess[162] = t283 * (t100 + t250);
    hess[163] = t283 * (-t107 + t108 + t265);
    hess[164] = -t283 * (-t114 + t276 + t7 * t9);
    hess[165] = t280;
    hess[166] = t281;
    hess[167] = t282;
    hess[168] =
        (1.0 / 4.0) * t278 * (t277 + 2) * (t35 - t64 * t7) / std::pow(d, 2);
}

} // namespace ipc::autogen
