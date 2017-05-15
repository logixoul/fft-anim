#include "precompiled.h"
#include "ciextra.h"
#include "fftwrap.h"
#include "labcolor.h"
#include <math.h> // bessel functions
#if 0

int sz=512;
Array2D<Vec2f> mainIn(sz, sz);
Array2D<float> overlay(sz, sz);
Array2D<float> overlay2(sz, sz);
typedef Array2D<Vec3f> Image;
typedef FFT_T<float> FFT;
ci::gl::GlslProg shader;
Image pic;
Image picOut;
int wscale=1;

float mouseX, mouseY;

template<class T> int mySign(T const& t) { if(t<0)return -1;if(t>0) return 1; return 0; }
struct Mover{
	Vec2i pos;
	Mover() { pos=Vec2i::zero(); }
	void move()
	{
		int random = ci::Rand::randInt(4);
		switch(random){
		case 0: pos.x+=1;break;
		case 1: pos.x-=1;break;
		case 2: pos.y+=1;break;
		case 3: pos.y-=1;break;
		}
		if(ci::Rand::randInt(16)==0) { pos.x -= mySign(pos.x); pos.y -= mySign(pos.y); }
		if(pos == Vec2i(0, 0) /*|| pos == Vec2i(-1,0) || pos == Vec2i(0, -1) || pos == Vec2i(-1,-1)*/) move();
	}
};
gl::Texture tex;
Mover mover;
bool evolve=false;
struct SApp : AppBasic {
	void setup()
	{
		createConsole();
		shader = gl::GlslProg(loadFile("shader.vs"), loadFile("shader.fs"));
		tex = gl::Texture(sz, sz);
		for(int x = 0; x < sz; x++) {
			for(int y= 0; y < sz; y++) {
				float dx = abs(x - sz/2);
				float dy = abs(y - sz/2);
				mainIn(x, y) = Vec2f(max(0.0f, 2 * (sz >> 3) - dx - dy), 0);
			}
		}
		setWindowSize(2*wscale*sz+3*10, wscale*sz+10);
		auto surface = Surface32f(loadImage("test.png"));
		pic = Image(surface.getWidth(), surface.getHeight());
		picOut = Image(surface.getWidth(), surface.getHeight());
		forxy(pic) pic(p) = ((Vec4f&)surface.getPixel(p)).xyz();
	}
	bool keys[256];
	void keyDown(KeyEvent e) {
		keys[e.getChar()]=true;
		if(e.getChar()=='e')evolve=!evolve;
		if(e.getChar()=='v'&&e.isControlDown())
		{
			ImageSourceRef surfaceRef = ci::Clipboard::getImage();
			if(!surfaceRef)
				return;
			Surface32f surface(surfaceRef);
			
			forxy(pic)
			{
				pic(p) = ((Vec4f&)surface.getPixel(p)).xyz();
			}
		}
	}
	void keyUp(KeyEvent e) { keys[e.getChar()]=false; }
	void mouseDown(MouseEvent e)
	{
	}
	void drawMoverInOverlay()
	{
		auto physPos = mover.pos;
		if(physPos.x<0)physPos.x+=sz;
		if(physPos.y<0)physPos.y+=sz;
		overlay(physPos) = .1;
	}
	void uploadComplexImage(FFT::CArray const& arr)
	{
		static Image display(sz, sz);
		complexArrayToImage(arr, display);
		/*forxy(display) {
			auto& c = display(p);
			
			auto bw = Vec3f::one() * Vec3f(1.1, .13, -.27).dot(c);
			Lab bwLab = ToLab(bw);
			Lab cLab = ToLab(c);
			cLab.L = bwLab.L;
			
			Vec2f ab(cLab.a, cLab.b);
			ab.rotate(2*M_PI*getMousePos().x / (float)getWindowWidth());
			cLab.a = ab.x;
			cLab.b = ab.y;

			Vec3f cBack = ToRgb(cLab);
			c = cBack;
		}*/
		tex.bind();
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, display.w, display.h, GL_RGB, GL_FLOAT, display.data);
	}
	// http://www.docjar.com/docs/api/java/awt/Color.html#HSBtoRGB%28float,%20float,%20float%29
	Color HSBtoRGB(float hue, float saturation, float brightness) 
	{
			int r = 0, g = 0, b = 0;
			if (saturation == 0) {
				r = g = b = (int) (brightness * 255.0f + 0.5f);
			} else {
				float h = (hue - (float)floor(hue)) * 6.0f;
				float f = h - (float)floor(h);
				float p = brightness * (1.0f - saturation);
				float q = brightness * (1.0f - saturation * f);
				float t = brightness * (1.0f - (saturation * (1.0f - f)));
				switch ((int) h) {
				case 0:
					r = (int) (brightness * 255.0f + 0.5f);
					g = (int) (t * 255.0f + 0.5f);
					b = (int) (p * 255.0f + 0.5f);
					break;
				case 1:
					r = (int) (q * 255.0f + 0.5f);
					g = (int) (brightness * 255.0f + 0.5f);
					b = (int) (p * 255.0f + 0.5f);
					break;
				case 2:
					r = (int) (p * 255.0f + 0.5f);
					g = (int) (brightness * 255.0f + 0.5f);
					b = (int) (t * 255.0f + 0.5f);
					break;
				case 3:
					r = (int) (p * 255.0f + 0.5f);
					g = (int) (q * 255.0f + 0.5f);
					b = (int) (brightness * 255.0f + 0.5f);
					break;
				case 4:
					r = (int) (t * 255.0f + 0.5f);
					g = (int) (p * 255.0f + 0.5f);
					b = (int) (brightness * 255.0f + 0.5f);
					break;
				case 5:
					r = (int) (brightness * 255.0f + 0.5f);
					g = (int) (p * 255.0f + 0.5f);
					b = (int) (q * 255.0f + 0.5f);
					break;
				}
			}
			return Color(r,g,b)/255.0f;
	}

	// from the applet
	Color complex_to_color(float re, float im, float R) {
		const float lmin=0.0, lmax=1.0, bmin=0.0, bmax=1.0;
		float hue, sat, brightness;
		float r;
		r = re * re + im * im;
		r = (float) sqrt(r);
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
	void complexArrayToImage(FFT::CArray const& in, Image& out)
	{
		const float sliderMax = getWindowHeight();
		const float sliderValue = getWindowHeight()/2;//getMousePos().y;
		float R = (float)tan((double)(sliderMax-sliderValue)/sliderMax*M_PI/2.0);

		forxy(in) {
			out(p) = (Vec3f&)complex_to_color(in(p).x, in(p).y, R);
			//out(p) *= 2;
		}
	}

	static void normalizeFft(FFT::CArray& arr)
	{
		float norm = 1.0 / sqrt((float)arr.area);
		for(int i=0; i < arr.area; i++) arr(i) *= norm;
	}
	static float Rh(float f) { return f/(f+1); }
	struct Gaussian
	{
		float valmul,expmul;
		Gaussian(float sigma) {
			float c = sigma;
			expmul=1.0f/(2*c*c);
			valmul=1.0/(sigma*sqrt(2*M_PI));
			valmul=1.0;
		}
		float eval(float f) { return valmul*exp(-f*f*expmul); }
	};
	
	void draw()
	{
		mouseX = getMousePos().x/(float)AppBasic::get()->getWindowWidth();
		mouseY = getMousePos().y/(float)AppBasic::get()->getWindowHeight();
		Gaussian gaussian(1/mouseY);
		gl::clear(Color(0, 0, 0));
		gl::setMatricesWindow(getWindowSize());
		glScalef(wscale, wscale, 1);
		if(evolve){
			for(int i = 0; i < 1; i++) {
				mover.move();
				drawMoverInOverlay();
			}
			forxy(overlay) {
				overlay(p) *= .9;
			}
			forxy(overlay) {
				overlay2(p) += overlay(p);
				overlay2(p) *= .9;
			}
		}
		tex.setMagFilter(GL_NEAREST);
		tex.setWrap(GL_REPEAT, GL_REPEAT);
		static FFT::CArray fftOut(sz, sz);
		static FFT::CArray fftBack(sz, sz);
		//cout<<"power: " << mouseX*10.0f << endl;
		static FFT::CArray sdKernel(sz, sz);
		static FFT::CArray fdKernel(sz, sz);
		forxy(sdKernel) {
			auto p2=p;if(p2.x>sz/2)p2.x-=sz;if(p2.y>sz/2)p2.y-=sz;
			float r = max(mouseX,2.0f/getWindowWidth());
			if(abs(getMousePos().x)<2) r = 1; // fix div by zero
			r*=100;
			//sdKernel(p).x=4*(1.0-smoothstep(r, r+1, p2.length()));
			//sdKernel(p).x-=3*(1.0-smoothstep(r-1, r, p2.length()));
			sdKernel(p).x=1.0-gaussian.eval(p2.length());
		}
		auto kernelInvSum = 1.0/(std::accumulate(sdKernel.begin(), sdKernel.end(), Vec2f::zero()).x);
		forxy(sdKernel) { sdKernel(p) *= kernelInvSum; }
		static auto kernelPlan = FFT::c2c(sdKernel, fdKernel, FFTW_FORWARD, FFTW_MEASURE); kernelPlan.execute();
		//cout<<"power "<<mouseX*10<<endl;
		for(int channel = 0; channel < 3; channel++)
		{
			auto coef=.99;
			forxy(pic) {
				auto val=pic(p).ptr()[channel];
				//val/=1-coef*val;
				mainIn(p) = Vec2f(0, val);
			}
			static auto plan1 = FFT::c2c(mainIn, fftOut, FFTW_FORWARD, FFTW_MEASURE); plan1.execute();
			normalizeFft(fftOut);
			forxy(fftOut) {
				fftOut(p) *= max(0.0, 1.0 - overlay2(p) * 100);
				auto& f = fftOut(p);
				//float f_lenSquared = f.lengthSquared();
				//float f2_len = pow(f_lenSquared, mouseX*10.0f-.5f);
				//f *= f2_len;//f *= (f2_len/f_len);
				//(complex<float>&)f*=(complex<float>&)fdKernel(p);
				auto p2=p;if(p2.x>sz/2)p2.x-=sz;if(p2.y>sz/2)p2.y-=sz;
				/*if(!Area(-1,-1,1,1).contains(p2))
					f *= 1.0-gaussian.eval(p2.length());*/
				float arg=p2.length();
				//f *= lerp((1.0/pow(arg+1, mouseX*10)), 1.0, smoothstep(0, 100, arg));
				//f*=(1.0/pow(arg+1, mouseX*10));
				const float e=2.71828182846;
				f*=(1.0/pow(log(arg+e), mouseX*10));
				//if(!Area(-1,-1,1,1).contains(p2))
				//	f*=smoothstep(1,20, arg);
				//f.rotate((int)(mouseX*10)*atan2(f.y,f.x));
			}
			static auto plan2 = FFT::c2c(fftOut, fftBack, FFTW_BACKWARD, FFTW_MEASURE); plan2.execute();
			normalizeFft(fftBack);
			forxy(picOut) {
				auto val= fftBack(p).length();
				val=max(val,0.0f);
				//val/=val+1;
				//val=pow(val, mouseY*10);
				picOut(p).ptr()[channel] = val;
			}
		}
		auto logfunc=[](float result, Vec3f const& v){return log(v.dot(Vec3f::one()/3));};
		auto logsum=[&](Image const& img) { return accumulate(img.begin(), img.end(), 0.0f, logfunc); };
		float scale = 1.0/exp(logsum(picOut));
		cout<<"scale"<<scale<<endl;
		forxy(picOut) picOut(p) *= scale;
		forxy(picOut) picOut(p) = apply(picOut(p), Rh);
		for(int channel=0; channel<3; channel++)
		{
			auto pred=[&](Vec3f const& a, Vec3f const& b){return a.ptr()[channel]<b.ptr()[channel];};
			float minval=0*std::min_element(picOut.begin(), picOut.end(), pred)->ptr()[channel];
			float maxval=std::max_element(picOut.begin(), picOut.end(), pred)->ptr()[channel];
			float mulval=1.0/(maxval-minval);
			forxy(picOut) { auto& c=picOut(p).ptr()[channel]; c-=minval; c*=mulval; }
		}
		float power=exp(mouseY*10);
		forxy(picOut) picOut(p) = apply(picOut(p), [&](float f){return pow(f, power);});
		//forxy(picOut) {picOut(p)+=pic(p)*.01;}
		
			
		/*shader.bind();
		shader.uniform("mouse", (Vec2f)getMousePos());
		shader.uniform("windowSize", (Vec2f)getWindowSize());*/
		//uploadComplexImage(fftBack);
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, picOut.w, picOut.h, GL_RGB, GL_FLOAT, picOut.data);
		gl::draw(tex, Rectf(10, 10, 10+sz, 10+sz));
		uploadComplexImage(fftOut);
		glMatrixMode(GL_TEXTURE); glPushMatrix(); glTranslatef(.5,.5,0.0); glMatrixMode(GL_MODELVIEW);
		gl::draw(tex, Rectf(sz+20, 10, sz+20+sz, 10+sz));
		glMatrixMode(GL_TEXTURE); glPopMatrix(); glMatrixMode(GL_MODELVIEW);
		::Sleep(150);
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