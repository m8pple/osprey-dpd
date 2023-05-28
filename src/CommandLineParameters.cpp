#include "StdAfx.h"
#include "SimDefs.h"
#include "DebugAssert.h"
#include "CommandLineParameters.h"

#include "xxCommandObject.h"

#include "IIntegrationEngine.h"

bool CommandLineParameters::sm_initialised = false;
std::vector<std::unique_ptr<xxCommandObject>> CommandLineParameters::sm_extra_commands;

static const char *sg_usage =
R"(
Usage: dpd [parameters] <runId1> <runId2>? <runId3>?

Executes runId1, then runId2, ...
If no runId is specified, then it will prompt the user to type one in.

The runId is the suffix of the dmpci file. For example, given "dmpci.water"
the runId would be "water". Note that the dmpci file should be in the current
working directory, e.g.:

# Assume there are two files in local directory
$ ls dmpci.*
dmpci.water dmpci.oil

# Run simulation in dmpci.water
$ dpd water

# Run simulation in dmpci.oil
$ dpd oil

# Run dmpci.water then dmpci.oil
$ dpd water oil

Parameters:

--add-command STRING : Interprets the given command as a Command, and inserts in the right place in sequence.

    STRING should be valid as a Command line in a dmpci file, though it is
    optional whether the "Command" token is included. The time must always be included. 
    Commands added in this way are not saved into the restart state, so if the simulation is restarted they
    must also be specified during the restart commands.

    Examples:

    # Run simulation in dmpci.water, setting time step size to 0.1 at time 4
    $ dpd --add-command "SetTimeStepSize 4 0.1" water

    # Equivalent to above, with optional Command prefix
    $ dpd --add-command "Command SetTimeStepSize 4 0.1" water

--set-engine EngineName : Selects a specific execution engine

--list-engines : Lists all available engines, then quits with error without starting simulation.

--wrap-engine-with-ref-diff : (Only useful for engine development) Wraps the current engine with
    a diff engine that compares bead state against the refernce engine output after each step.
    If no engine is set then it is a no-op.

)";

void CommandLineParameters::Initialise(int &argc, char **&argv, std::function<void(const std::string &)> on_error)
{
    auto consume_args = [&](unsigned n)
    {
        DEBUG_ASSERT( n < argc );
        for(unsigned i=1; i<argc-n; i++){
            argv[i] = argv[i+n];
        }
        argc -= n;
    };

    if(sm_initialised){
        return on_error("CommandLineParameters::Initialise already called.");
    }

    sm_initialised=true;

    while(argc > 1){
        std::string cmd=argv[1];
        if(cmd.size() < 2){
            break;
        }
        if(cmd.substr(0,2) != "--"){
            break;
        }

        if(cmd == "--add-command"){
            if(argc < 3){
                return on_error("Missing argument to --add-command");
            }
            std::string command_string = argv[2] + std::string("\n");
            std::stringstream src(command_string);
            std::string command_name;
            long command_time;
            src >> command_name >> command_time;
            if(!src.good()){
                return on_error("Coudn't get command name and time from '"+command_string+"'");
            }

            std::unique_ptr<xxCommandObject> cmd(new xxCommandObject(command_name, command_time));
            src >> *cmd;
            if(!src.good())
            {
                return on_error("Error reading Command data for " + command_name + " from '"+command_string+"'");
            }
            if(!cmd->IsCommandValid())
            {
                return on_error("Command not value for " + command_name + " from '"+command_string+"'");
            }

            std::cout << *cmd;

            sm_extra_commands.push_back(std::move(cmd));

            consume_args(2);
        }else if(cmd=="--list-engines"){
            IIntegrationEngineFactory::ListEngines(std::cout);
            exit(1);
        }else if(cmd=="--set-engine"){
            if(argc < 3){
                return on_error("Missing argument to --set-engine");
            }

            std::string engine_name = argv[2];
            
            auto inst = IIntegrationEngineFactory::CreateEngineInstanceByName(engine_name);
            if(!inst){
                return on_error("Couldn't find engine called "+engine_name);
            }
            IIntegrationEngine::SetGlobalEngine(inst );

            consume_args(2);
        }else if(cmd=="--wrap-engine-with-ref-diff"){
            IIntegrationEngine::WrapGlobalEngineWithRefDiff();
            consume_args(1);
        }else if(cmd=="--help"){
            consume_args(1);
            return on_error(sg_usage);
        }else{
            return on_error("Didn't understand command line flag '"+cmd+"'");
        }
    }
}