#include "precompiled.h"
#include "ciextra.h"
#include "fftwrap.h"
#include "labcolor.h"
#include "hdrwrite.h"
#if 1
int sz=256;
Array2D<Vec2f> mainIn(sz, sz);
Array2D<Vec2f> mainInBackup(sz, sz);
Array2D<float> overlay(sz, sz);
typedef Array2D<Vec3f> Image;
typedef FFT_T<float> FFT;
ci::gl::GlslProg shader;
gl::Texture mask,mask2;//,dest;
int wscale=2;
float exposure;

float mouseX, mouseY;
/*struct xyzrgb_struct {
	ci::Matrix33f rgb2xyz;
	ci::Matrix33f xyz2rgb;
	xyzrgb_struct()
	{
		rgb2xyz = (1.0 / 0.17697) *
			mat3(
			Vec3f(0.49, 0.17697, 0.00),
			Vec3f(0.31, 0.81240, 0.01),
			Vec3f(0.20, 0.01063, 0.99)
			);
		rgb2xyz.inverted
		xyz2rgb = inverse(rgb2xyz);
	}
} xyzrgb;*/

template<class T> int mySign(T const& t) { if(t<0)return -1;if(t>0) return 1; return 0; }
struct Mover{
	Vec2f pos;
	Vec2i center;
	int iterations;
	float power;
	Mover(Vec2i center, int iterations, int power) {
		pos=Vec2f::zero(); this->center = center; this->iterations = iterations;
		this->power = power;
	}
	void move()
	{
		int random = ci::Rand::randInt(4);
		float s = ci::randFloat()*.01;///(float)iterations;
		switch(random){
		case 0: pos.x+=s;break;
		case 1: pos.x-=s;break;
		case 2: pos.y+=s;break;
		case 3: pos.y-=s;break;
		}
		pos -= center;
		pos*=.999f;
		pos += center;
		//if((Vec2i)pos == Vec2i(0, 0) /*|| pos == Vec2i(-1,0) || pos == Vec2i(0, -1) || pos == Vec2i(-1,-1)*/) move();
	}
};
gl::Texture tex;
gl::Fbo fbo;
vector<Mover> movers;
bool evolve=true;
struct SApp : AppBasic {
	void setup()
	{
		createConsole();
		shader = gl::GlslProg(loadFile("shader.vs"), loadFile("shader.fs"));
		gl::Texture::Format fmt1;
		fmt1.setInternalFormat(GL_RGBA32F);
		gl::Fbo::Format rgba32f;
		rgba32f.setColorInternalFormat(GL_RGBA32F);
		gl::Texture::Format srgb;
		srgb.setInternalFormat(GL_SRGB8_ALPHA8_EXT);
		tex = gl::Texture(sz, sz, fmt1);
		for(int x = 0; x < sz; x++) {
			for(int y= 0; y < sz; y++) {
				float dx = abs(x - sz/2);
				float dy = abs(y - sz/2);
				mainInBackup(x, y) = Vec2f(max(0.0f, 2 * (sz >> 3) - dx - dy), 0)*0.0f;
			}
		}
		mainIn = mainInBackup.clone();
		setWindowSize(2*wscale*sz+3*10, wscale*sz+10);
		mask = gl::Texture(loadImage("mask.png"), srgb);
		mask2 = gl::Texture(loadImage("mask2.png"), srgb);
		movers.push_back(Mover(Vec2i(0, 0), 1, 1));
		movers.push_back(Mover(Vec2i(0, 0), 10, 10));
		movers.push_back(Mover(Vec2i(0, 0), 100, 100));
		movers.push_back(Mover(Vec2i(0, 0), 1000, 1000));
		fbo=gl::Fbo(wscale*sz, wscale*sz, rgba32f);
		ci::Rand::randSeed(11);
	}
	bool keys[256];
	void keyDown(KeyEvent e) {
		keys[e.getChar()]=true;
		if(e.getChar()=='e')evolve=!evolve;
		if(e.getChar()=='s'&&e.isControlDown())
			save();
	}
	void save()
	{
		static int saveIndex=0;
		Image readback(fbo.getSize());
		fbo.bindTexture();
		glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_FLOAT, readback.data);
		auto f = fopen(("hdr output/"+boost::lexical_cast<string>(saveIndex)+".hdr").c_str(),"wb");
		rgbe_header_info header;
		header.valid=RGBE_VALID_EXPOSURE;
		header.exposure=exposure;
		RGBE_WriteHeader(f,readback.w,readback.h,&header);
		RGBE_WritePixels(f,(float*)readback.data,readback.area);
		fclose(f);
		saveIndex++;
	}
	void keyUp(KeyEvent e) { keys[e.getChar()]=false; }
	void mouseDown(MouseEvent e)
	{
	}
	// http://www.docjar.com/docs/api/java/awt/Color.html#HSBtoRGB%28float,%20float,%20float%29
	Color HSBtoRGB(float hue, float saturation, float brightness) 
	{
			float r = 0, g = 0, b = 0;
			if (saturation == 0) {
				r = g = b = brightness;
			} else {
				float h = (hue - (float)floor(hue)) * 6.0f;
				float f = h - (float)floor(h);
				float p = brightness * (1.0f - saturation);
				float q = brightness * (1.0f - saturation * f);
				float t = brightness * (1.0f - (saturation * (1.0f - f)));
				switch ((int) h) {
				case 0:
					r =  (brightness);
					g =  (t);
					b =  (p);
					break;
				case 1:
					r =  (q);
					g =  (brightness);
					b =  (p);
					break;
				case 2:
					r =  (p);
					g =  (brightness);
					b =  (t);
					break;
				case 3:
					r =  (p);
					g =  (q);
					b =  (brightness);
					break;
				case 4:
					r =  (t);
					g =  (p);
					b =  (brightness);
					break;
				case 5:
					r =  (brightness);
					g =  (p);
					b =  (q);
					break;
				}
			}
			return Color(r,g,b);
	}
	Color HSBtoRGB2(float hue, float saturation) 
	{
			float r = 0, g = 0, b = 0;
			if (saturation == 0) {
				return Color(1.0f, 1.0f, 1.0f);
			} else {
				float h = (hue - (float)floor(hue)) * 6.0f;
				float f = h - (float)floor(h);
				float p = 1.0f - saturation;
				float q = 1.0f - saturation * f;
				float t = 1.0f - saturation * (1.0f - f);
				switch ((int) h) {
				case 0:
					return Color(1.0f, t, p);
					break;
				case 1:
					return Color(q, 1.0f, p);
					break;
				case 2:
					return Color(p, 1.0f, t);
					break;
				case 3:
					return Color(p, q, 1.0f);
					break;
				case 4:
					return Color(t, p, 1.0f);
					break;
				case 5:
					return Color(1.0f, p, q);
					break;
				}
			}
	}
	void uploadComplexImage(FFT::CArray const& arr)
	{
		static Image display(sz, sz);

		complexArrayToImage(arr, display);
		tex.bind();
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, display.w, display.h, GL_RGB, GL_FLOAT, display.data);
	}
	// from the applet
	Color complex_to_color(float re, float im, float R) {
		const float lmin=0.0, lmax=1.0, bmin=0.0, bmax=1.0;
		float hue, sat, brightness;
		float r;
		r = re * re + im * im;
		//r = (float) sqrt(r);
		float li = std::min<float>(1.0f, (std::max<float>(0.0f,
				(float)atan(r/R) * 2.0f/(float)M_PI)*(bmax-bmin) + bmin)) *
				(lmax-lmin)+lmin;
		if (li<0.5f) {
		  sat = 1.0f;
		  brightness = 2.f*li;
		}
		else {
		  sat = 2.0f-2.0f*li;
		  brightness = 1.0f;
		}
		//brightness = ((float)atan(r/R) * 2f/(float)PI);
		float phi = (float)M_PI+(float)atan2(im,re);
		double scale = (double)4.0/20.0;
		hue = (float)( phi+scale*sin(3*phi) ) / (float)(2*M_PI);  // Correct for the "banding" in HSB color where RGB are wider than the others
		return HSBtoRGB(hue, sat, brightness);
	}
	Color complex_to_color2(float re, float im) {
		float hue, sat, brightness;
		float r;
		r = re * re + im * im;
		r = (float) sqrt(r);
		//brightness = ((float)atan(r/R) * 2f/(float)PI);
		float phi = (float)M_PI+(float)atan2(im,re);
		//double scale = (double)4.0/20.0;
		hue = (float)( phi /*+scale*sin(3*phi)*/ ) / (float)(2*M_PI);  // Correct for the "banding" in HSB color where RGB are wider than the others
		//return (Color&)(hsvToRGB(Vec3f(hue, 1.0f, 1.0f)) * atan(brightness/R));
		
		auto c = HSBtoRGB2(hue, .99f) /* *r */;
		return c * exp(r * 6.0f);
		
		/*Color rgb=HSBtoRGB(hue, .95f, 1.0f);
		Vec3f vec(rgb.r,rgb.g,rgb.b);
		vec.safeNormalize();
		vec*=r;
		return Color(vec.x, vec.y, vec.z);*/
	}
	void complexArrayToImage(FFT::CArray const& in, Image& out)
	{
		/*const float sliderMax = getWindowHeight();
		const float sliderValue = getMousePos().y;
		float R = (float)tan((double)(sliderMax-sliderValue)/sliderMax*M_PI/2.0);
		*/
		forxy(in) {
			out(p) = (Vec3f&)complex_to_color2(in(p).x, in(p).y);
		}
	}
	void drawMoverInFftOut(FFT::CArray& arr, Mover const& mover)
	{
		auto moverPos2 = mover.pos * mover.power;
		auto physPos = moverPos2;
		if(physPos.x<0)physPos.x+=sz;
		if(physPos.y<0)physPos.y+=sz;
		if(!Area(0, 0, sz, sz).contains(physPos))
			return;
		(complex<float>&)arr(physPos) += .00001f*(complex<float>&)ci::Rand::randVec2f() / pow(1.0f+moverPos2.length(), 1.0f);
	}
	void draw()
	{
		auto blendc=Vec3f(182.0, 160.0, 194.0)/255.0;
		blendc = Vec3f::one() / (Vec3f::one() - (blendc-Vec3f::one()*.5)*2.0);
		glClampColorARB(GL_CLAMP_FRAGMENT_COLOR_ARB, GL_FALSE);
		glClampColorARB(GL_CLAMP_VERTEX_COLOR_ARB, GL_FALSE);
		mouseX = getMousePos().x/(float)AppBasic::get()->getWindowWidth();
		mouseY = getMousePos().y/(float)AppBasic::get()->getWindowHeight();
		gl::clear(Color(0, 0, 0));
		gl::setMatricesWindow(getWindowSize());
		glScalef(wscale, wscale, 1);
		//tex.setMagFilter(GL_NEAREST);
		mask.setWrap(GL_CLAMP, GL_CLAMP);
		static FFT::CArray originalFftOut(sz, sz);
		static FFT::CArray fftOut(sz, sz);
		static FFT::CArray fftBack(sz, sz);
		static auto plan1 = FFT::c2c(mainIn, fftOut, FFTW_FORWARD, FFTW_MEASURE); plan1.execute();
		for(int i=0; i < fftOut.area; i++) fftOut(i) /= sqrt((float)fftOut.area); // normalize
		static bool first = true;
		if(first) {
			originalFftOut=fftOut.clone();
			first=false;
		}
		if(evolve){
			foreach(auto& mover, movers) {
				for(int i = 0; i < mover.iterations; i++) {
					mover.move();
					drawMoverInFftOut(fftOut, mover);
				}
			}
			forxy(fftOut) {
				fftOut(p) = lerp(fftOut(p), originalFftOut(p), .0000001);
			}
		}
		forxy(fftOut) {
			fftOut(p) *= max(0.0, 1.0 - overlay(p) * 100);
			/*auto& f = fftOut(p);
			auto f2 = (Vec2f&)pow((complex<float>&)f, mouseX*10.0f);
			f = f2.safeNormalized() * f.length();*/
		}
		static auto plan2 = FFT::c2c(fftOut, fftBack, FFTW_BACKWARD, FFTW_MEASURE); plan2.execute();
		for(int i=0; i < fftBack.area; i++) fftBack(i) /= sqrt((float)fftBack.area); // normalize
		float magSum=0.0f; forxy(fftBack) magSum+=fftBack(p).length();
		if(magSum!=0.0){
			float magmul=1.0/(magSum/(float)fftBack.area);
			forxy(fftBack) fftBack(p)*=magmul;
			magmul = 1.0 / magmul;
			uploadComplexImage(fftBack);
			forxy(fftBack) fftBack(p)*=magmul;
		}
		std::copy(fftBack.begin(), fftBack.end(), mainIn.begin());
		shader.bind();
		shader.uniform("mouse", (Vec2f)getMousePos());
		shader.uniform("windowSize", (Vec2f)getWindowSize());
		shader.uniform("tex", 0); tex.bind(0);
		shader.uniform("mask", 1); mask.bind(1); 
		shader.uniform("mask2", 2); mask2.bind(2);
		shader.uniform("time", (float)getElapsedSeconds());
		//glColor3f(exposure,exposure,exposure);
		tex.setWrap(GL_CLAMP, GL_CLAMP);
		fbo.bindFramebuffer();
		//gl::draw(tex, Rectf(10, 10, 10+sz, 10+sz));
		gl::pushMatrices(); glLoadIdentity(); glPushAttrib(GL_ALL_ATTRIB_BITS); glViewport(0, 0, fbo.getWidth(), fbo.getHeight());
		gl::setMatricesWindow(fbo.getSize());
		gl::draw(tex, fbo.getBounds());
		glPopAttrib(); gl::popMatrices();
		fbo.unbindFramebuffer();
		shader.unbind();
		gl::draw(fbo.getTexture(), Rectf(10, 10, 10+sz, 10+sz));
		
		/*uploadComplexImage(fftOut);
		glColor3f(1e2, 1e2, 1e2);
		glMatrixMode(GL_TEXTURE); glPushMatrix(); glTranslatef(.5,.5,0.0); glMatrixMode(GL_MODELVIEW);
		tex.setWrap(GL_REPEAT, GL_REPEAT);
		gl::draw(tex, Rectf(sz+20, 10, sz+20+sz, 10+sz));
		glMatrixMode(GL_TEXTURE); glPopMatrix(); glMatrixMode(GL_MODELVIEW);
		glColor3f(1,1,1);*/
	}
};

//CINDER_APP_BASIC(SApp, RendererGl)
int WINAPI WinMain(HINSTANCE hInstance,HINSTANCE hPrevInstance,LPSTR lpCmdLine,int nCmdShow) {	
	try{
		cinder::app::AppBasic::prepareLaunch();														
		cinder::app::AppBasic *app = new SApp;														
		cinder::app::Renderer *ren = new RendererGl;													
		cinder::app::AppBasic::executeLaunch( app, ren, "SApp" );										
		cinder::app::AppBasic::cleanupLaunch();														
		return 0;																					
	}catch(ci::gl::GlslProgCompileExc const& e) {
		cout << "caught: " << endl << e.what() << endl;
		//int dummy;cin>>dummy;
		system("pause");

	}
}
#endif