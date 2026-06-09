#include "agent.h"

static size_t WriteCallback(void *contents, size_t size, size_t nmemb, void *userp)
{
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

Agent::Agent(Grid grid) {
    this->grid = grid;
    stepSize = 0.05 * grid.getEdges()[0]->getCurrentLength();
}

json Agent::statAgent() {
    return runCURL("http://127.0.0.1:8000/stat_model");
}

json Agent::runCURL(string url) {
    curl = curl_easy_init();
    if(!curl) {
        cout << "Error: cURL failed to initialize\n";
        exit(1);
    }
    readBuffer = "";

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
    struct curl_slist *hs=NULL;
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, hs);

    res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);

    auto out = json::parse(readBuffer);
    return out;
}

json Agent::runCURL(json &inputs) {
    return runCURL(inputs, "http://127.0.0.1:8000/run_model");
}

json Agent::runCURL(json &inputs, string url) {
    curl = curl_easy_init();
    if(!curl) {
        cout << "Error: cURL failed to initialize\n";
        exit(1);
    }
    readBuffer = "";

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
    struct curl_slist *hs=NULL;
    hs = curl_slist_append(hs, "Content-Type: application/json");
    curl_easy_setopt(curl, CURLOPT_POST, 1L);
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, hs);

    auto payload = inputs.dump();
    curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, payload.length());
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, payload.c_str());
    res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);

    auto out = json::parse(readBuffer);
    return out;
}

json Agent::runAgent(json &inputs) {
    auto results = runCURL(inputs);
    return results;
}
