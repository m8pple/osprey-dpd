#ifndef FastOstream_h
#define FastOstream_h

#include <cmath>
#include <vector>
#include "DebugAssert.hpp"
#include "xxState.h"



/*
An extremely minimal formatter which only knows how to format:
- Integers
- Floating-point numbers in fixed-point format with dst.precision() decimal digits
- chars
- strings
*/
struct FastOstream
{
	std::vector<char> buffer;
	char *curr=0;
	
	double scale;
	int digits;

	xxState *context;
	std::ostream &dst;

	FastOstream(xxState *_context, std::ostream &_dst)
		: context(_context)
		, dst(_dst)
	{
		buffer.resize(4096);
		curr=&buffer[0];

		digits=dst.precision();
		assert(digits>0);
		// scale=std::exp10(digits); OSX libm doesn't have this ? Anyway, linker errors...
		scale=1;
		for(int i=0; i<digits; i++){
			scale *= 10;
		}
	}

	FastOstream(const FastOstream &) = delete;
	FastOstream &operator=(const FastOstream &)=delete;
	
	~FastOstream()
	{
		flush();
	}

	FastOstream &operator<<(long x)
	{
        if(x < 0){
            x = -x;
            *curr++='-';
        }

		if(x<10){ // Handles 0 and other single digit cases. Common for bead types
			*curr++='0'+x;
		}else{
			static_assert(sizeof(long) <= 8);
			char scratch[21]; // A long can't be more than this
			char * begin=scratch + sizeof(scratch);
            char *end=begin;

			while(x>0){
				*--begin = '0' + (x%10);
				x = x/10; 
			}

            DEBUG_ASSERT(begin < end);
            memcpy(curr, begin, end-begin);
            curr += end-begin;
		}

		return *this;
	}

	FastOstream &operator<<(char ch)
	{
		*curr++ = ch;
		if(ch=='\n'){
			flush();
		}
		return *this;
	}

	FastOstream &operator<<(double x)
	{
		int c=std::fpclassify(x);
		switch(c){
		case FP_INFINITE:
        case FP_NAN:
			flush();
			dst << x;
            return *this;
        case FP_NORMAL:
            break;
        case FP_SUBNORMAL:
        case FP_ZERO:
            *curr++ = '0';
			return *this;
        default:
            DEBUG_ASSERT(0);
            if(context){    
    			context->FatalTrace("Invalid floating point classification.");
            }else{
                xxState::FatalTraceGlobal("Invalid floating-point classification.");
            }
		}

		if(x < 0){
			*curr++ = '-';
			x = -x;
		}

		if( x > 1000000.0 ){
			flush();
			dst << x;
		}

		static_assert(sizeof(long) <= 8);
		char scratch[32]; // A fixed-point double in this range can't be more tahn this
        char *begin=scratch+sizeof(scratch);
        char *end=begin;

		// Scale to fixed point with target digit count
		int64_t xs = std::round(x * scale);

		// unconditionally output the fractional digits
		for(int i=0; i<digits; i++){
			*--begin = '0' + xs % 10;
			xs = xs / 10;
		}

		// Decimal point
		*--begin = '.';

		// Unconditionally output first decimal
		*--begin = '0' + xs % 10;
		xs = xs / 10;

		// If more than on decimal then keep going
		while(xs > 0){
			*--begin = '0' + xs % 10;
			xs = xs / 10;
		}

        memcpy(curr, begin, end-begin);
        curr += end-begin;

		return *this;
	}



	void flush()
	{
		if(curr==&buffer.back()){
			return;
		}

		assert(&buffer[0] <= curr);
		assert(curr < &buffer.back());
		if(curr >= &buffer.back()){
            if(context){
    			context->FatalTrace("Internal error : a fast formatted line used too many characters.");
            }else{
                xxState::FatalTraceGlobal("Internal error : a fast formatted line used too many characters.");
            }
		}
		dst.write(&buffer[0], curr - &buffer[0]);
		curr=&buffer[0];
	}

};

#endif
