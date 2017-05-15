#version 130
const float pi = 3.14159265358979323846264;

struct Lab { float L, a, b; };
struct Lch { float L, c, h; };
struct Xyz { float X, Y, Z; };

float refX =  95.047;
float refY = 100.000;
float refZ = 108.883;

Lab ToLab(Xyz xyz) {
	float x = xyz.X / refX;
    float y = xyz.Y / refY;
    float z = xyz.Z / refZ;
    
    if ( x > 0.008856 ) x = pow(x, 1.0/3.0);
    else                    x = ( 7.787 * x ) + ( 16.0 / 116.0 );
    if ( y > 0.008856 ) y = pow(y, 1.0/3.0);
    else                    y = ( 7.787 * y ) + ( 16.0 / 116.0 );
    if ( z > 0.008856 ) z = pow(z, 1.0/3.0);
    else                    z = ( 7.787 * z ) + ( 16.0 / 116.0 );
    
	return Lab(
        ( 116.0 * y ) - 16.0,
		500.0 * ( x - y ),
		200.0 * ( y - z )
	    );
}

Xyz ToXyz(Lab lab) {
	float y = ( lab.L + 16.0 ) / 116.0;
	float x = lab.a / 500.0 + y;
	float z = y - lab.b / 200.0;

	if ( pow(y, 3.0) > 0.008856 ) y = pow(y, 3.0);
	else                      y = ( y - 16.0 / 116.0 ) / 7.787;
	if ( pow(x, 3.0) > 0.008856 ) x = pow(x, 3.0);
	else                      x = ( x - 16.0 / 116.0 ) / 7.787;
	if ( pow(z, 3.0) > 0.008856 ) z = pow(z, 3.0);
	else                      z = ( z - 16.0 / 116.0 ) / 7.787;

	return Xyz(refX * x, refY * y, refZ * z);
}

vec3 ToRgb(Xyz xyz) {
    float x = xyz.X / 100.0;
    float y = xyz.Y / 100.0;
    float z = xyz.Z / 100.0;
    float var_R = x * 3.2406 + y * -1.5372 + z * -0.4986;
    float var_G = x * -0.9689 + y * 1.8758 + z * 0.0415;
    float var_B = x * 0.0557 + y * -0.2040 + z * 1.0570;

    if (var_R > 0.0031308) var_R = 1.055 * pow(var_R, 1.0 / 2.4) - 0.055;
    else var_R = 12.92 * var_R;
    if (var_G > 0.0031308) var_G = 1.055 * pow(var_G, 1.0 / 2.4) - 0.055;
    else var_G = 12.92 * var_G;
    if (var_B > 0.0031308) var_B = 1.055 * pow(var_B, 1.0 / 2.4) - 0.055;
    else var_B = 12.92 * var_B;

    return vec3(var_R, var_G, var_B);
}

Xyz ToXyz(vec3 rgb) {
    float var_R = rgb.r;
    float var_G = rgb.g;
    float var_B = rgb.b;

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

Lch ToLch(Lab lab)
{
    float var_H = atan(lab.b, lab.a);

    if ( var_H > 0.0 ) var_H = ( var_H / pi ) * 180.0;
    else             var_H = 360.0 - ( abs( var_H ) / pi ) * 180.0;

    return Lch(lab.L, length(vec2(lab.a, lab.b)), var_H);
}

Lab ToLab(Lch lch)
{
    return Lab(
        lch.L,
        cos(radians(lch.h)) * lch.c,
        sin(radians(lch.h)) * lch.c
        );
}

Lab ToLab(vec3 rgb) { return ToLab(ToXyz(rgb)); }
vec3 ToRgb(Lab lab) { return ToRgb(ToXyz(lab)); }
Lch ToLch(vec3 rgb) { return ToLch(ToLab(rgb)); }
vec3 ToRgb(Lch lch) { return ToRgb(ToLab(lch)); }

//////////

mat2 rotationMatrix(float angle) {
	return mat2( cos( angle ), -sin( angle ),
		sin( angle ),  cos( angle ));
}

#define BlendColorDodgef(base, blend) 	((blend == 1.0) ? blend : min(base / (1.0 - blend), 1.0))
#define BlendColorBurnf(base, blend) 	((blend == 0.0) ? blend : max((1.0 - ((1.0 - base) / blend)), 0.0))
#define BlendVividLightf(base, blend) 	((blend < 0.5) ? BlendColorBurnf(base, (2.0 * blend)) : BlendColorDodgef(base, (2.0 * (blend - 0.5))))
#define Blend(base, blend, funcf) 		vec3(funcf(base.r, blend.r), funcf(base.g, blend.g), funcf(base.b, blend.b))
#define BlendVividLight(base, blend) 	vec3(BlendVividLightf(base.r, blend.r), BlendVividLightf(base.g, blend.g), BlendVividLightf(base.b, blend.b))
vec3 RGBToHSL(vec3 color)
{
	vec3 hsl; // init to 0 to avoid warnings ? (and reverse if + remove first part)
	
	float fmin = min(min(color.r, color.g), color.b);    //Min. value of RGB
	float fmax = max(max(color.r, color.g), color.b);    //Max. value of RGB
	float delta = fmax - fmin;             //Delta RGB value

	hsl.z = (fmax + fmin) / 2.0; // Luminance

	if (delta == 0.0)		//This is a gray, no chroma...
	{
		hsl.x = 0.0;	// Hue
		hsl.y = 0.0;	// Saturation
	}
	else                                    //Chromatic data...
	{
		if (hsl.z < 0.5)
			hsl.y = delta / (fmax + fmin); // Saturation
		else
			hsl.y = delta / (2.0 - fmax - fmin); // Saturation
		
		float deltaR = (((fmax - color.r) / 6.0) + (delta / 2.0)) / delta;
		float deltaG = (((fmax - color.g) / 6.0) + (delta / 2.0)) / delta;
		float deltaB = (((fmax - color.b) / 6.0) + (delta / 2.0)) / delta;

		if (color.r == fmax )
			hsl.x = deltaB - deltaG; // Hue
		else if (color.g == fmax)
			hsl.x = (1.0 / 3.0) + deltaR - deltaB; // Hue
		else if (color.b == fmax)
			hsl.x = (2.0 / 3.0) + deltaG - deltaR; // Hue

		if (hsl.x < 0.0)
			hsl.x += 1.0; // Hue
		else if (hsl.x > 1.0)
			hsl.x -= 1.0; // Hue
	}

	return hsl;
}

float HueToRGB(float f1, float f2, float hue)
{
	if (hue < 0.0)
		hue += 1.0;
	else if (hue > 1.0)
		hue -= 1.0;
	float res;
	if ((6.0 * hue) < 1.0)
		res = f1 + (f2 - f1) * 6.0 * hue;
	else if ((2.0 * hue) < 1.0)
		res = f2;
	else if ((3.0 * hue) < 2.0)
		res = f1 + (f2 - f1) * ((2.0 / 3.0) - hue) * 6.0;
	else
		res = f1;
	return res;
}

vec3 HSLToRGB(vec3 hsl)
{
	vec3 rgb;
	
	if (hsl.y == 0.0)
		rgb = vec3(hsl.z); // Luminance
	else
	{
		float f2;
		
		if (hsl.z < 0.5)
			f2 = hsl.z * (1.0 + hsl.y);
		else
			f2 = (hsl.z + hsl.y) - (hsl.y * hsl.z);
			
		float f1 = 2.0 * hsl.z - f2;
		
		rgb.r = HueToRGB(f1, f2, hsl.x + (1.0/3.0));
		rgb.g = HueToRGB(f1, f2, hsl.x);
		rgb.b= HueToRGB(f1, f2, hsl.x - (1.0/3.0));
	}
	
	return rgb;
}

float clamp01(float f) { return clamp(f, 0.0, 1.0); }

//imagemagick
vec3 ConvertHCLToRGB(float hue, float chroma, float luma)
{
  float b,c,g,h,m,r,x,z;

  /*
    Convert HCL to RGB colorspace.
  */
  h=6.0*hue;
  c=chroma;
  x=c*(1.0-abs(mod(h,2.0)-1.0));
  r=0.0;
  g=0.0;
  b=0.0;
  if ((0.0 <= h) && (h < 1.0))
    {
      r=c;
      g=x;
    }
  else
    if ((1.0 <= h) && (h < 2.0))
      {
        r=x;
        g=c;
      }
    else
      if ((2.0 <= h) && (h < 3.0))
        {
          g=c;
          b=x;
        }
      else
        if ((3.0 <= h) && (h < 4.0))
          {
            g=x;
            b=c;
          }
        else
          if ((4.0 <= h) && (h < 5.0))
            {
              r=x;
              b=c;
            }
          else
            if ((5.0 <= h) && (h < 6.0))
              {
                r=c;
                b=x;
              }
  m=luma-(0.298839*r+0.586811*g+0.114350*b);
  /*
    Choose saturation strategy to clip it into the RGB cube; hue and luma are
    preserved and chroma may be changed.
  */
  z=1.0;
  if (m < 0.0)
    {
      z=luma/(luma-m);
      m=0.0;
    }
  else
    if (m+c > 1.0)
      {
        z=(1.0-luma)/(m+c-luma);
        m=1.0-z*c;
      }
		return vec3(clamp01(z*r+m), clamp01(z*g+m), clamp01(z*b+m));
}

vec3 ConvertRGBToHCL(vec3 rgb)
{
  float b,c,g,h,max,r;
  float red=rgb.r;
  float green=rgb.g;
  float blue=rgb.b;

  /*
    Convert RGB to HCL colorspace.
  */
  r=rgb.r;
  g=rgb.g;
  b=rgb.b;
  max=max(r,max(g,b));
  c=max-min(r,min(g,b));
  h=0.0;
  if (c == 0)
    h=0.0;
  else
    if (red == max)
      h=mod(6.0+(g-b)/c,6.0);
    else
      if (green == max)
        h=((b-r)/c)+2.0;
      else
        if (blue == max)
          h=((r-g)/c)+4.0;
  float hue=(h/6.0);
  float chroma=c;
  float luma=0.298839*r+0.586811*g+0.114350*b;
  return vec3(hue, chroma, luma);
}

// Luminosity Blend mode creates the result color by combining the hue and saturation of the base color with the luminance of the blend color.
vec3 BlendLuminosity(vec3 base, vec3 blend)
{
	vec3 baseHCL = ConvertRGBToHCL(base);
	return ConvertHCLToRGB(baseHCL.x, baseHCL.y, ConvertRGBToHCL(blend).z);
}

varying vec2 tc;
uniform sampler2D tex;
uniform sampler2D mask;
uniform sampler2D mask2;
uniform vec2 mouse;
uniform vec2 windowSize;
uniform float time;

vec3 Filmic(vec3 c)
{
	vec3 x = max(vec3(0.0),c-0.004);
	return (x*(6.2*x+.5))/(x*(6.2*x+1.7)+0.06);
}

vec3 mod289(vec3 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 mod289(vec4 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 permute(vec4 x) {
     return mod289(((x*34.0)+1.0)*x);
}

vec4 taylorInvSqrt(vec4 r)
{
  return 1.79284291400159 - 0.85373472095314 * r;
}

float snoise(vec3 v)
  { 
  const vec2  C = vec2(1.0/6.0, 1.0/3.0) ;
  const vec4  D = vec4(0.0, 0.5, 1.0, 2.0);

// First corner
  vec3 i  = floor(v + dot(v, C.yyy) );
  vec3 x0 =   v - i + dot(i, C.xxx) ;

// Other corners
  vec3 g = step(x0.yzx, x0.xyz);
  vec3 l = 1.0 - g;
  vec3 i1 = min( g.xyz, l.zxy );
  vec3 i2 = max( g.xyz, l.zxy );

  //   x0 = x0 - 0.0 + 0.0 * C.xxx;
  //   x1 = x0 - i1  + 1.0 * C.xxx;
  //   x2 = x0 - i2  + 2.0 * C.xxx;
  //   x3 = x0 - 1.0 + 3.0 * C.xxx;
  vec3 x1 = x0 - i1 + C.xxx;
  vec3 x2 = x0 - i2 + C.yyy; // 2.0*C.x = 1/3 = C.y
  vec3 x3 = x0 - D.yyy;      // -1.0+3.0*C.x = -0.5 = -D.y

// Permutations
  i = mod289(i); 
  vec4 p = permute( permute( permute( 
             i.z + vec4(0.0, i1.z, i2.z, 1.0 ))
           + i.y + vec4(0.0, i1.y, i2.y, 1.0 )) 
           + i.x + vec4(0.0, i1.x, i2.x, 1.0 ));

// Gradients: 7x7 points over a square, mapped onto an octahedron.
// The ring size 17*17 = 289 is close to a multiple of 49 (49*6 = 294)
  float n_ = 0.142857142857; // 1.0/7.0
  vec3  ns = n_ * D.wyz - D.xzx;

  vec4 j = p - 49.0 * floor(p * ns.z * ns.z);  //  mod(p,7*7)

  vec4 x_ = floor(j * ns.z);
  vec4 y_ = floor(j - 7.0 * x_ );    // mod(j,N)

  vec4 x = x_ *ns.x + ns.yyyy;
  vec4 y = y_ *ns.x + ns.yyyy;
  vec4 h = 1.0 - abs(x) - abs(y);

  vec4 b0 = vec4( x.xy, y.xy );
  vec4 b1 = vec4( x.zw, y.zw );

  //vec4 s0 = vec4(lessThan(b0,0.0))*2.0 - 1.0;
  //vec4 s1 = vec4(lessThan(b1,0.0))*2.0 - 1.0;
  vec4 s0 = floor(b0)*2.0 + 1.0;
  vec4 s1 = floor(b1)*2.0 + 1.0;
  vec4 sh = -step(h, vec4(0.0));

  vec4 a0 = b0.xzyw + s0.xzyw*sh.xxyy ;
  vec4 a1 = b1.xzyw + s1.xzyw*sh.zzww ;

  vec3 p0 = vec3(a0.xy,h.x);
  vec3 p1 = vec3(a0.zw,h.y);
  vec3 p2 = vec3(a1.xy,h.z);
  vec3 p3 = vec3(a1.zw,h.w);

//Normalise gradients
  vec4 norm = taylorInvSqrt(vec4(dot(p0,p0), dot(p1,p1), dot(p2, p2), dot(p3,p3)));
  p0 *= norm.x;
  p1 *= norm.y;
  p2 *= norm.z;
  p3 *= norm.w;

// Mix final noise value
  vec4 m = max(0.6 - vec4(dot(x0,x0), dot(x1,x1), dot(x2,x2), dot(x3,x3)), 0.0);
  m = m * m;
  return 42.0 * dot( m*m, vec4( dot(p0,x0), dot(p1,x1), 
                                dot(p2,x2), dot(p3,x3) ) );
  }
float lengthSquared(vec2 v) { return dot(v, v); }

// begin new version of xyz code

// from cinder
mat3 inverseMat(mat3 m, float epsilon)
{
	mat3 inv;

	// Compute the adjoint.
	inv[0][0] = m[1][1]*m[2][2] - m[1][2]*m[2][1];
	inv[0][1] = m[0][2]*m[2][1] - m[0][1]*m[2][2];
	inv[0][2] = m[0][1]*m[1][2] - m[0][2]*m[1][1];
	inv[1][0] = m[1][2]*m[2][0] - m[1][0]*m[2][2];
	inv[1][1] = m[0][0]*m[2][2] - m[0][2]*m[2][0];
	inv[1][2] = m[0][2]*m[1][0] - m[0][0]*m[1][2];
	inv[2][0] = m[1][0]*m[2][1] - m[1][1]*m[2][0];
	inv[2][1] = m[0][1]*m[2][0] - m[0][0]*m[2][1];
	inv[2][2] = m[0][0]*m[1][1] - m[0][1]*m[1][0];

	float det = m[0][0]*inv[0][0] + m[0][1]*inv[1][0] + m[0][2]*inv[2][0];

	if( abs( det ) > epsilon ) {
		float invDet = 1.0/det;
		inv *= invDet;
		/*inv[0][0] *= invDet;
		inv[0][1] *= invDet;
		inv[0][2] *= invDet;
		inv[1][0] *= invDet;
		inv[1][1] *= invDet;
		inv[1][2] *= invDet;
		inv[2][0] *= invDet;
		inv[2][1] *= invDet;
		inv[2][2] *= invDet;*/
	}

	return inv;
}

//uniform mat3 rgb2xyz;
//uniform mat3 xyz2rgb;
mat3 rgb2xyz = (1.0 / 0.17697) *
	mat3(
	vec3(0.49, 0.17697, 0.00),
	vec3(0.31, 0.81240, 0.01),
	vec3(0.20, 0.01063, 0.99)
	);
// each arg is a column in the matrix
vec3 Xyz2Rgb(vec3 xyz) {
// first row of matrix is coefs for first output
	mat3 xyz2rgb = inverseMat(rgb2xyz, 0.0000001);
	return xyz2rgb*xyz;
}

vec3 Rgb2Xyz(vec3 rgb) {
    return rgb2xyz*rgb;
}

// end new version of xyz code

// begin YUV

mat3 rgb2yuv = mat3(
	vec3(0.299, -0.14713, 0.615),
	vec3(0.587, -0.28886, -0.51499),
	vec3(0.114, 0.436, -0.10001)
	);
mat3 yuv2rgb = inverseMat(rgb2yuv, 0.000001);

vec3 Rgb2Yuv(vec3 rgb) {
	return rgb2yuv*rgb;
}
vec3 Yuv2Rgb(vec3 yuv) {
	return yuv2rgb*yuv;
}

// end YUV

void main()
{
	vec2 tc2 = tc-vec2(.5);
	vec3 _out = vec3(0.0);
	float exponent=1.0;
	for(int i = 0; i < 80; i++)
	{
		exponent *= .9;
		tc2*=0.99;
		vec3 c = texture2D(tex, tc2+vec2(.5)).xyz;
		//if(c!=vec3(0.0))
			//c=normalize(c)*(1.0+exp(length(c)*8.0)/1000.0);
		
		c*=exponent;
		if(i!=0)
			c *= 0.01;
		_out += c;
	}
	float exposure=exp(-15+10.0*mouse.x/windowSize.x);
	_out*=exposure;

	vec3 xyz=Rgb2Yuv(_out);
	xyz.x /= xyz.x + 1.0;
	_out = Yuv2Rgb(xyz);

	//_out/=_out+vec3(1.0);
	_out=pow(_out, vec3(1.0/2.2));
	//_out = Filmic(_out);

	
	gl_FragColor.rgb = _out;
	gl_FragColor.a = 1.0;
}