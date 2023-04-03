#ifndef ISimEngine_h
#define ISimEngine_h

#include <string>
#include <memory>

#include <unordered_map>

class ISimBox;

class ISimEngineCapabilities
{
public:
    virtual ~ISimEngineCapabilities()
    {}

    virtual std::string Name() const=0;

    virtual bool IsParallel() const=0;

    virtual std::string CanSupport(const ISimBox *box) const =0;
};


class ISimEngine
    : public ISimEngineCapabilities
{
public:
    virtual ~ISimEngine()
    {}

    //! /param modified If true the some other process may have modified the box since the last step. Otherwise it is untouched
    virtual void Run(ISimBox *box, bool modified, unsigned num_steps) =0;


    static std::shared_ptr<ISimEngine> GetGlobalEngine()
    {
        return sim_engine_backing();
    }

    static void SetGlobalEngine(std::shared_ptr<ISimEngine> engine)
    {
        sim_engine_backing() = engine;
    }
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

private:
    static std::unordered_map<std::string,std::shared_ptr<ISimEngineFactory>> &GetFactories()
    {
        static std::unordered_map<std::string,std::shared_ptr<ISimEngineFactory>> factories;
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
        std::shared_ptr<ISimEngine> CreateEngineInstance()
        { return std::make_shared<TImpl>(); }

        std::string Name() const override
        { return TCapabilities::Name(); }

        bool IsParallel() const override
        { return TCapabilities::IsParallel(); }

        std::string CanSupport(const ISimBox *box) const override
        { return TCapabilities::CanSupport(box); }
    };

public:
    static bool Register()
    {
        ISimEngineFactory::RegisterEngine(std::make_shared<Factory>());
        return true;
    }
};

#endif
