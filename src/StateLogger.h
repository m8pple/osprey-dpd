#ifndef StateLogger_hpp
#define StateLogger_hpp

#include <string>
#include <array>
#include <functional>

#if defined(OSPREY_DPD_ENABLE_STATE_LOGGING)
constexpr bool StateLogger_CompileTimeEnable = 1;
#else
constexpr bool StateLogger_CompileTimeEnable = 0;
#endif

class StateLogger
{
public:
    // Definitive false if disabled at compile time, so all logging code
    // optimised out.
    static bool IsEnabled()
    {
        return StateLogger_CompileTimeEnable && m_RunTimeEnabled;
    }

    static void Enable(std::function<void(const std::string &)> line_sink);
    static void Disable();


    static void SetPrefix(const std::string &_prefix);
    static void BeginStep(unsigned step);


    // These helper functions are in the header so that if logging is
    // disabled at run-time they will be completely optimised out.

    static void Log(const char *facet, int key_dim, const int *key, int val_dim, const double *value)
    {
        if(StateLogger_CompileTimeEnable){
            LogImpl(facet, key_dim, key, val_dim, value);
        }
    }

    static void Log(const char *facet, double value)
    {
        if(StateLogger_CompileTimeEnable){
            LogImpl(facet, 0, 0, 1, &value);
        }
    }

    // Log a scalar property for b1
    static void LogBead(const char *facet, unsigned b1,  double value)
    {
        if(StateLogger_CompileTimeEnable){
            int bb=b1;
            Log(facet, 1, &bb, 1, &value);
        }
    }

    // Log a vector property for b1
    template<class T>
    static void LogBead(const char *facet, unsigned b1, const T *value)
    {
        if(StateLogger_CompileTimeEnable){
            int bb=b1;
            double v[3]={value[0],value[1],value[2]};
            Log(facet, 1, &bb, 3, v);
        }
    }

    // Log a shared scalar property just for (b1,b2)
    static void LogBeadPair(const char *facet, unsigned b1, unsigned b2, double value)
    {
        if(StateLogger_CompileTimeEnable){
            int bb[2]={(int)b1,(int)b2};
            Log(facet, 2, bb, 1, &value);
        }
    }

    // Log a shared vector property just for (b1,b2)
    template<class T>
    static void LogBeadPair(const char *facet, unsigned b1, unsigned b2, const T *value)
    {
        if(StateLogger_CompileTimeEnable){
            int bb[2]={(int)b1,(int)b2};
            double pos[3]={value[0], value[1], value[2]};
            Log(facet, 2, bb, 3, pos);
        }
    }

    // Log a shared vector property just for (b1,b2)
    static void LogBeadPair(const char *facet, unsigned b1, unsigned b2, const std::array<double,3> &values)
    {
        if(StateLogger_CompileTimeEnable){
            int bb[2]={(int)b1,(int)b2};
            Log(facet, 2, bb, 3, &values[0]);
        }
    }

    // Log a shared scalar property for (b1,b2), and the same property for (b2,b1) 
    static void LogBeadPairRefl(const char *facet, unsigned b1, unsigned b2, double value)
    {
        if(StateLogger_CompileTimeEnable){
            int bb[2]={(int)b1,(int)b2};
            Log(facet, 2, bb, 1, &value);
            std::swap(bb[0], bb[1]);
            Log(facet, 2, bb, 1, &value);
        }
    }

    // Log a vector property for (b1,b2), and the negated property for (b2,b1)
    template<class T>
    static void LogBeadPairRefl(const char *facet, unsigned b1, unsigned b2, const T *value)
    {
        if(StateLogger_CompileTimeEnable){
            int bb[2]={(int)b1,(int)b2};
            double pos[3]={value[0], value[1], value[2]};
            Log(facet, 2, bb, 3, pos);
            std::swap(bb[0], bb[1]);
            double neg[3]={-value[0], -value[1], -value[2]};
            Log(facet, 2, bb, 3, neg);
        }
    }

    static void LogBeadPairRefl(const char *facet, unsigned b1, unsigned b2, const std::array<double,3> &value)
    {
        if(StateLogger_CompileTimeEnable){
            LogBeadPairRefl(facet, b1, b2, &value[0]);
        }
    }


    // Log a shared scalar property for (b1,b2,b3)
    static void LogBeadTriple(const char *facet, unsigned b1, unsigned b2, unsigned b3, double value)
    {
        if(StateLogger_CompileTimeEnable){
            int bb[3]={(int)b1,(int)b2, (int)b3};
            Log(facet, 3, bb, 1, &value);
        }
    }

    // Log a shared vector property for (b1,b2,b3)
    template<class T>
    static void LogBeadTriple(const char *facet, unsigned b1, unsigned b2, unsigned b3, const T *value)
    {
        if(StateLogger_CompileTimeEnable){
            int bb[3]={(int)b1,(int)b2,(int)b3};
            double v[3]={value[0],value[1], value[2]};
            Log(facet, 3, bb, 3, value);
        }
    }

    // Log a shared vector property for (b1,b2,b3)
    static void LogBeadTriple(const char *facet, unsigned b1, unsigned b2, unsigned b3, const std::array<double,3> &values)
    {
        if(StateLogger_CompileTimeEnable){
            int bb[3]={(int)b1,(int)b2,(int)b3};
            Log(facet, 3, bb, 3, &values[0]);
        }
    }

private:
    static void LogImpl(const char *facet, int key_dim, const int *key, int val_dim, const double *value);

    static bool m_RunTimeEnabled;
};

#endif
