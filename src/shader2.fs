uniform sampler2D tex; uniform sampler2D blurred; varying vec2 tc;
uniform vec2 mouse;
uniform vec2 blurredScale;
uniform vec2 tex_tsize;
uniform sampler2D env;

const float pi = 3.14159265358979323846264;

vec3 Filmic(vec3 c)
{
	vec3 x = max(vec3(0.0),c-0.004);
	return (x*(6.2*x+.5))/(x*(6.2*x+1.7)+0.06);
}

vec3 w = vec3(.22, .71, .07);

void clipChannel(inout vec3 v, int channel)
{
	if(v[channel] > 1.0)
	{
		vec3 equalLgray = vec3(dot(v, w));
		vec3 p = equalLgray - v; // coefs of linear equation parametric representation
		// | result[channel] = 1
		// | result(s) = v + s * p
		// | result[channel] = v[channel] + s * p[channel] = 1
		// | s = (1-v[channel]) / p[channel]
		float s = (1.0-v[channel]) / p[channel];
		vec3 result = v + s * p;
		v = result;
	}
}

vec2 tex_tsize2 = tex_tsize * .5;
float g_ = 5.0;

vec3 get_(vec2 offset)
{
	vec3 a = texture2D(tex, tc + tex_tsize2 * offset).rgb;
	return a / (a + vec3(1.0));
}

vec3 getSharpenedMainTex(out vec3 normal)
{
	vec3 a00 = get_(vec2(-1.0, -1.0));
	vec3 a01 = get_(vec2(-1.0, +0.0));
	vec3 a02 = get_(vec2(-1.0, +1.0));
	vec3 a10 = get_(vec2(+0.0, -1.0));
	vec3 a11 = get_(vec2(+0.0, +0.0));
	vec3 a12 = get_(vec2(+0.0, +1.0));
	vec3 a20 = get_(vec2(+1.0, -1.0));
	vec3 a21 = get_(vec2(+1.0, +0.0));
	vec3 a22 = get_(vec2(+1.0, +1.0));
	float f = .3;
	vec3 b = max(vec3(0.0), mix(a11, 9.0 * a11 - (a00 + a01 + a02 + a10 + a12 + a20 + a21 + a22), f));
	b = b / (vec3(1.0) - b);

	vec3 v1 = vec3(1.0 / 3.0);
	float dx = dot(a01, v1) - dot(a21, v1);
	float dy = dot(a10, v1) - dot(a12, v1);
	float b2 = sqrt(dot(texture2D(blurred, tc * blurredScale).rgb, v1));
	dx /= b2;
	dy /= b2;
	normal = normalize(vec3(dx, dy, 1.0));

	return b;
}

void main() {
	
	vec3 N;
	vec3 _out = getSharpenedMainTex(N);
	vec3 L3 = normalize(vec3(0.0, 0.0, 1.2));
	vec3 I = normalize(vec3(tc.x, tc.y, 0.0) - vec3(.5, .5, .3));
	//N = vec3(0.0, 0.0, 1.0);
	vec3 R = reflect(I, N);
	//R = vec3(R.x, -R.z, R.y);
	_out += texture2D(env, vec2(atan(R.x, R.z) / (2.0 * pi) + .5 + .5, asin(-R.y) / pi + .5)).rgb;
	//if(dot(_out, vec3(1.0/3.0)) > .01)
	//	_out += /*dot(_out, vec3(1.0/3.0)) * */ vec3(pow(max(0.0, dot(R, L3)), 40.0)) * 10.0;
	vec3 b = texture2D(blurred, tc * blurredScale).rgb * exp(-10.0+10.0*mouse.y);
	//vec3 b2 = pow(b, vec3(2.0));
	//b = (b2 / dot(b2, w)) * dot(b, w);
	//_out += b;

	
	if(mouse.y < 0.0) {
	float L2 = dot(_out, w);
	_out = mix(_out, vec3(L2), L2 / (L2 + 1.0)); }
	
	//_out += b;

	float L = dot(_out, w);
	_out /= L;
	_out *= L / (L + 1.0);

	if(mouse.y >= 0.0) {
		clipChannel(_out, 0);
		clipChannel(_out, 1);
		clipChannel(_out, 2);
	}
	
	//_out = vec3(1.0) - (vec3(1.0)-b) * (vec3(1.0)-_out);
	_out = mix(_out, vec3(1.0), b);
	//_out += b;

	_out = pow(_out, vec3(1.0/2.2));
	//_out *= 3.0;
	//_out = Filmic(_out);
	gl_FragColor = vec4(_out, 1.0);
}