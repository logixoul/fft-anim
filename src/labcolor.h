#pragma once
#include "precompiled.h"

struct Lab { float L, a, b; Lab(float L, float a, float b):L(L),a(a),b(b){} };
struct Lch { float L, c, h; Lch(float L, float c, float h):L(L),c(c),h(h){} };
struct Xyz { float X, Y, Z; Xyz(float X, float Y, float Z):X(X),Y(Y),Z(Z){} };
Lab ToLab(Xyz const& xyz);
Xyz ToXyz(Lab const& lab);
Vec3f ToRgb(Xyz const& xyz);
Xyz ToXyz(Vec3f const& rgb);

inline Lab ToLab(Vec3f const& rgb) { return ToLab(ToXyz(rgb)); }
inline Vec3f ToRgb(Lab const& lab) { return ToRgb(ToXyz(lab)); }