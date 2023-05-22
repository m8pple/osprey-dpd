#ifndef ISimEngine_h
#define ISimEngine_h

#include <string>
#include <memory>
#include <iostream>

#include <map>

class ISimBox;

class ISimEngineCapabilities
{
public:
    virtual ~ISimEngineCapabilities()
    {}

    virtual std::string Name() const=0;

    virtual bool IsParallel() const=0;

    enum SupportStatus
    {
        Supported = 0,

        // Problems that may clear if the sim box is stepped at least one step with a different engine
        TransientProblemStep,

        // Problems that may clear if a flag is turned off. Effectively
        // they need a command to change the simbox state.
        TransientProblemFlags,

        // Problems that are inherent to the sim box or compilation
        PermanentProblem
    };

    struct support_result
    {
        support_result(const support_result &) = default;
        support_result(support_result &&) = default;
        support_result &operator=(const support_result &) = default;
        support_result &operator=(support_result &&) = default;

        support_result(SupportStatus _status)
            : status(_status)
        {
            assert(_status==Supported);
        }

        support_result(SupportStatus _status, const std::string &_reason)
            : status(_status)
            , reason(_reason)
        {}

        SupportStatus status;
        std::string reason; 
    };

    /* Check whether the engine can run this box.
        If the box is supported then the status==Supported and reason is empty.
        Otherwise the status tells you whether it is a transient problem
        or a permanent fault, and the reason gives at least one reason why.

        Transient faults may clear 
    */
    support_result CanSupport(const ISimBox *box) const;

    /*
    These are standard assumptions required by all SimEngine instances. They
    are properties of compilation and the sim box, so if they are not met
    they will never be met.
    
    Compile time (PermanentProblem):
    - SimDimensions == 3
    - EnableStressTensorSphere is not enabled
    - EnableParallelSimBox is not enabled (i.e. not an MPI build)
    - UseDPDBeadRadii is not defined (all bead radii == 1)
    Run time SimBox (PermanentProblem):
    - lambda==0.5
    - Box origin is at zero
    - Box lengths are integers
    - Cell widths are integers
    Run time flags (TransientProblem):
    - Active commands to run in current time-step
    - 
    */
    static support_result StandardSupportAssumptions(const ISimBox *box);
protected:
    // Used by engines to add to the standard constraints
    virtual support_result CanSupport_ExtraConstraints(const ISimBox *box) const
    {
        return {Supported, {}};
    }

};


class ISimEngine
    : public ISimEngineCapabilities
{
public:
    virtual ~ISimEngine()
    {}

    struct run_result
    {
        run_result(SupportStatus _status, std::string _reason, unsigned _steps)
            : status(_status)
            , reason(_reason)
            , completed_steps(_steps)
        {
            assert( (_status != Supported) == _reason.size()!=0 );
        }

        run_result(unsigned _steps)
            : status(Supported)
            , completed_steps(_steps)
        {}

        run_result(support_result res)
            : status(res.status)
            , reason(res.reason)
            , completed_steps(0)
        {
            assert( (status != Supported) == reason.size()!=0 );
        }

        SupportStatus status;
        std::string reason; 
        unsigned completed_steps;
    };

    /*! Attempts to run the box for the given number of time-steps
    \retval Returns the number of steps actually run.
    \param modified If true the some other process may have modified the box since the last step. Otherwise it is untouched
    
    The assumption is that most of the time the engine will be running the box,
    with only statistical monitoring steps in-between. This is the case where
    the modified flag is false - the _only_ things that might have changed are
    the positions, momentum, and forces of the beads.

    If anything else has changed (new beads, beads change type, polymers changed, ...)
    the modified flag should be set to allow the engine to rebuild any data structures
    and/or check it can still support it.

    Engines may run less steps than requested for various reasons, such as:
    - The box is not supported. Though ideally the called would check this.
    - Some internal assumptions were broken, e.g.:
        - beads are moving too fast for the algorithm
        - the local cell density got too high.
    In such cases the engine will guarantee that the box state is preserved,
    either by not modifying it, rolling it back, or some other process. When
    there are temporary violations and engine might try to go as far as possible,
    then stop at the problem point. Another engine could then try to take
    over. These are reflected in the run result code:
    - Success : requested number of steps was completed
    - TransientFault : less than the expected number was completed, but if it is
        stepped in another engine it might be able to carry on
    - PermanentFault : This engine cannot support the box for some reason 
    */
    virtual run_result Run(ISimBox *box, bool modified, unsigned start_sim_time, unsigned num_steps) =0;


    static std::shared_ptr<ISimEngine> GetGlobalEngine()
    {
        return sim_engine_backing();
    }

    static void SetGlobalEngine(std::shared_ptr<ISimEngine> engine)
    {
        sim_engine_backing() = engine;
    }

    /*! Wraps the current engine (if any) with a diff engine that checks the results against
        the reference engine after every time step.
    */
    static void WrapGlobalEngineWithRefDiff();
private:
    static std::shared_ptr<ISimEngine> &sim_engine_backing()
    {
        static std::shared_ptr<ISimEngine> backing;
        return backing;
    }
};

class ISimEngineFactory
    : public ISimEngineCapabilities
{
public:
    virtual std::shared_ptr<ISimEngine> CreateEngineInstance() =0;

    static std::shared_ptr<ISimEngine> CreateEngineInstanceByName(std::string name)
    {
        auto &factories=GetFactories();
        auto it=factories.find(name);
        if(it==factories.end()){
            return {};
        }

        return it->second->CreateEngineInstance();
    }

    static bool RegisterEngine(std::shared_ptr<ISimEngineFactory> factory)
    {
        auto &factories=GetFactories();
        if(factories.find(factory->Name()) != factories.end()){
            fprintf(stderr, "Error : SimEngineFactories for %s added twice\n", factory->Name().c_str());
            return false;
        }

        factories[factory->Name()] = factory;

        return true;
    }

    static bool IsKnownEngine(const std::string &name)
    {
        auto &factories=GetFactories();
        return factories.find(name) != factories.end();
    }

    static void ListEngines(std::ostream &dst)
    {
        for(const auto &kv : GetFactories()){
            dst<<kv.first<<"\n";
        }
    }

private:
    static std::map<std::string,std::shared_ptr<ISimEngineFactory>> &GetFactories()
    {
        static std::map<std::string,std::shared_ptr<ISimEngineFactory>> factories;
        return factories;
    }
};


template<class TImpl, class TCapabilities=TImpl>
class SimEngineBase
{
private:
    class Factory
        : public ISimEngineFactory
        , public TCapabilities
    {
        std::shared_ptr<ISimEngine> CreateEngineInstance() override
        { return std::make_shared<TImpl>(); }

        std::string Name() const override
        { return TCapabilities::Name(); }

        bool IsParallel() const override
        { return TCapabilities::IsParallel(); }

        ISimEngineCapabilities::support_result CanSupport_ExtraConstraints(const ISimBox *box) const override
        { return TCapabilities::CanSupport_ExtraConstraints(box); }
    };

public:
    static bool Register()
    {
        ISimEngineFactory::RegisterEngine(std::make_shared<Factory>());
        return true;
    }
};

#endif
