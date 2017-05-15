#include "precompiled.h"
#include "labcolor.h"

float refX =  95.047;
float refY = 100.000;
float refZ = 108.883;

Lab ToLab(Xyz const& xyz) {
	float x = xyz.X / refX;
    float y = xyz.Y / refY;
    float z = xyz.Z / refZ;
    
    if ( x > 0.008856 ) x = pow(x, 1.0f/3.0f);
    else                    x = ( 7.787 * x ) + ( 16.0 / 116.0 );
    if ( y > 0.008856 ) y = pow(y, 1.0f/3.0f);
    else                    y = ( 7.787 * y ) + ( 16.0 / 116.0 );
    if ( z > 0.008856 ) z = pow(z, 1.0f/3.0f);
    else                    z = ( 7.787 * z ) + ( 16.0 / 116.0 );
    
	return Lab(
        ( 116.0 * y ) - 16.0,
		500.0 * ( x - y ),
		200.0 * ( y - z )
	    );
}

Xyz ToXyz(Lab const& lab) {
	float y = ( lab.L + 16.0 ) / 116.0;
	float x = lab.a / 500.0 + y;
	float z = y - lab.b / 200.0;

	if ( pow(y, 3.0f) > 0.008856 ) y = pow(y, 3.0f);
	else                      y = ( y - 16.0 / 116.0 ) / 7.787;
	if ( pow(x, 3.0f) > 0.008856 ) x = pow(x, 3.0f);
	else                      x = ( x - 16.0 / 116.0 ) / 7.787;
	if ( pow(z, 3.0f) > 0.008856 ) z = pow(z, 3.0f);
	else                      z = ( z - 16.0 / 116.0 ) / 7.787;

	return Xyz(refX * x, refY * y, refZ * z);
}

Vec3f ToRgb(Xyz const& xyz) {
    float x = xyz.X / 100.0;
    float y = xyz.Y / 100.0;
    float z = xyz.Z / 100.0;
    float var_R = x * 3.2406 + y * -1.5372 + z * -0.4986;
    float var_G = x * -0.9689 + y * 1.8758 + z * 0.0415;
    float var_B = x * 0.0557 + y * -0.2040 + z * 1.0570;

    if (var_R > 0.0031308) var_R = 1.055 * pow(var_R, 1.0f / 2.4f) - 0.055;
    else var_R = 12.92 * var_R;
    if (var_G > 0.0031308) var_G = 1.055 * pow(var_G, 1.0f / 2.4f) - 0.055;
    else var_G = 12.92 * var_G;
    if (var_B > 0.0031308) var_B = 1.055 * pow(var_B, 1.0f / 2.4f) - 0.055;
    else var_B = 12.92 * var_B;

    return Vec3f(var_R, var_G, var_B);
}

Xyz ToXyz(Vec3f const& rgb) {
    float var_R = rgb.x;
    float var_G = rgb.y;
    float var_B = rgb.z;

    if (var_R > 0.04045) var_R = pow((var_R + 0.055) / 1.055, 2.4);
    else var_R = var_R / 12.92;
    if (var_G > 0.04045) var_G = pow((var_G + 0.055) / 1.055, 2.4);
    else var_G = var_G / 12.92;
    if (var_B > 0.04045) var_B = pow((var_B + 0.055) / 1.055, 2.4);
    else var_B = var_B / 12.92;

    var_R = var_R * 100.0;
    var_G = var_G * 100.0;
    var_B = var_B * 100.0;

    return Xyz(
        var_R * 0.4124 + var_G * 0.3576 + var_B * 0.1805,
        var_R * 0.2126 + var_G * 0.7152 + var_B * 0.0722,
        var_R * 0.0193 + var_G * 0.1192 + var_B * 0.9505
        );
}
