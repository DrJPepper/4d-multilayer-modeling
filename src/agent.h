#pragma once

#include "util.h"
#include "grid.h"

class Agent {
    public:
        Agent(Grid);
        ~Agent() = default;
        json runAgent(json &inputs);
        json statAgent();

    private:
        json runCURL(json &input);
        json runCURL(json &input, string url);
        json runCURL(string url);
        bool onCPU;
        CURL *curl;
        CURLcode res;
        std::string readBuffer;
        num stepSize;
        Grid grid;
};
