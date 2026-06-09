#include <toml++/impl/parse_error.hpp>
#define DOCTEST_CONFIG_IMPLEMENT
//#define RUN_DOCTEST
#include "src/util.h"
#include "src/widget.h"
#include "src/simulation.h"

toml::table setts;
ColVector3d transFinal;
Matrix3d rotFinal;

int main(int argc, char *argv[]) {

#ifdef RUN_DOCTEST
    doctest::Context ctx;
    ctx.setOption("abort-after", 500);  // default - stop after 5 failed asserts
    ctx.applyCommandLine(argc, argv); // apply command line - argc / argv
    ctx.setOption("no-breaks", true); // override - don't break in the debugger
    int res = ctx.run();              // run test cases unless with --no-run
    if(ctx.shouldExit())              // query flags (and --exit) rely on this
        return res;
    exit(0);
#endif

    string inputFileName;
    // Parse arguments
    // Declare the supported options.
    po::options_description desc("Usage: ./agency_sim input-file [additional args]");
    desc.add_options()
        ("help,h", "Produce help message")
        ("dont-load-save,d", "Do not load save file of model even if present")
        ("input-file,f", po::value<string>(), "(Required) TOML configuration file location")
        ("import-from-json,i", "Import saved grid from a JSON file")
        ("scene-saving,j", "Generate JSON versions of scenes for export")
        ("headless,l", "Run in headless mode")
    ;
    po::positional_options_description p;
    p.add("input-file", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(desc).positional(p).allow_unregistered().run(), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        cout << desc << "\n";
        return 0;
    }
    if (vm.count("input-file")) {
        inputFileName = vm["input-file"].as<string>();
    } else {
        cerr << "ERROR: input-file argument is required\n";
        return 1;
    }

    std::ifstream inFile(inputFileName);
    if (!inFile.good()) {
        cerr << "ERROR: input file invalid or doesn't exist\n";
        return 1;
    }

    try {
        setts = toml::parse_file(inputFileName);
    } catch (const toml::parse_error& err) {
        cerr << "Parsing TOML file failed:\n" << err << "\n";
        return 1;
    }
    // TODO: Make it so these can also be put in the TOML file directly
    toml::table setts2;
    setts2.insert("import_from_json", static_cast<bool>(vm.count("import-from-json")));
    setts2.insert("headless", static_cast<bool>(vm.count("headless")));
    setts2.insert("scene_saving", static_cast<bool>(vm.count("scene-saving")));
    setts2.insert("dont_load_save", static_cast<bool>(vm.count("dont-load-save")));
    setts.insert("global", setts2);

    Simulation *sim;
    string t = setts["grid"]["type"].value_or("epi");
    sim = new Simulation(!t.compare("cubic"));
    if (!setts["global"]["headless"].value_or(false)) {
        // needed to ensure appropriate OpenGL context is created for VTK rendering.
        QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());

        QApplication a(argc, argv);
        string fn;
        if (setts["global"]["import_from_json"])
            fn = inputFileName;
        else
            fn = sim->getModel().modelFile;
        Widget *w = new Widget(0, sim, fn);
        w->show();
        return a.exec();
    } else {
        sim->runHeadless();
    }
    return 0;
}
