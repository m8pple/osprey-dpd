#include "StateLogger.h"
#include <string>
#include <sstream>
#include <functional>

bool StateLogger::m_RunTimeEnabled = false;

static unsigned StateLogger_Stepindex =  0;
static std::string StateLogger_Prefix;
static std::string StateLogger_PrefixPlusTime;
static std::function<void(const std::string &)> StateLogger_LineSink;
static std::mutex StateLogger_Mutex;


void StateLogger::Enable(std::function<void(const std::string &)> line_sink)
{
    StateLogger_LineSink = line_sink;
    if(line_sink){
        StateLogger::m_RunTimeEnabled=true;
    }else{
        StateLogger::m_RunTimeEnabled=false;
    }
}

void StateLogger::Disable()
{
    StateLogger_LineSink = {};
    StateLogger::m_RunTimeEnabled=false;
}


void StateLogger::SetPrefix(const std::string &_prefix)
{
    StateLogger_Prefix = _prefix;
    StateLogger_PrefixPlusTime=StateLogger_Prefix+std::to_string(StateLogger_Stepindex)+",";
}

void StateLogger::BeginStep(unsigned step)
{
    StateLogger_Stepindex = step;
    StateLogger_PrefixPlusTime=StateLogger_Prefix+std::to_string(StateLogger_Stepindex)+",";
}

void StateLogger::LogImpl(const char *facet, int key_dim, const int *key, int val_dim, const double *value)
{
    if(IsEnabled()){
        std::stringstream line;
        line<<StateLogger_PrefixPlusTime;
        line<<facet;
        for(int dim=0; dim<3; dim++){
            line<<',';
            if(dim<key_dim){
                line<<key[dim];
            }
        }
        for(int dim=0; dim<3; dim++){
            line<<',';
            if(dim<val_dim){
                line<<value[dim];
            }
        }
        {
            std::lock_guard<std::mutex> lk(StateLogger_Mutex);
            StateLogger_LineSink(line.str());
        }
    }
}
