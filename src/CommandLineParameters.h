#ifndef CommandLineParameters_h
#define CommandLineParameters_h

#include <vector>
#include <functional>
#include <memory>

class xxCommandObject;


class CommandLineParameters
{
public:
    static const std::vector<std::unique_ptr<xxCommandObject>> &GetExtraCommands()
    { return sm_extra_commands; }

    static void Initialise(int &argc, char **&argv, std::function<void(const std::string &)> on_error);
private:
    static bool sm_initialised;
    static std::vector<std::unique_ptr<xxCommandObject>> sm_extra_commands;
};

#endif