// bro.mesh / bro.rigging under the host's access gate (embed/embed_access.h;
// src/api/api.h, "Access"): every loader asks Read and every saver Write about
// the path it was handed, a refusal is an Error carrying the gate's sentence,
// code EACCES and the absolute path, nothing is read or written, and with the
// gate removed everything behaves as before. Run from test_mesh_api.cpp's main.

#include "eval/eval.h"
#include "embed/embed.h"
#include "embed/embed_access.h"

#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace {

namespace fs = std::filesystem;
namespace ev = bronze::embed;

#define GCHECK(cond, msg)                                                      \
    do {                                                                       \
        if (!(cond)) {                                                         \
            std::cerr << "CHECK failed: " #cond " (" << __FILE__ << ":"        \
                      << __LINE__ << "): " << (msg) << std::endl;              \
            std::exit(1);                                                      \
        }                                                                      \
    } while (0)

std::string utf8(const fs::path& p) {
    const std::u8string u = p.u8string();
    return std::string(u.begin(), u.end());
}

// The path as the gate is handed it.
std::string gatePath(const fs::path& p) {
    return utf8(fs::absolute(p).lexically_normal());
}

std::string jsStr(const fs::path& p) {
    std::string out = "'";
    for (char c : utf8(p)) {
        if (c == '\\' || c == '\'') out += '\\';
        out += c;
    }
    return out + "'";
}

std::string evalString(const std::string& src) {
    auto r = bronze::eval::evalScript(src);
    if (r.thrown) return "<thrown>";
    return ev::isString(r.value) ? ev::toUtf8(r.value) : "<not a string>";
}

// Runs `stmt`; "ok" or the refusal as one comparable string.
std::string syncJs(const std::string& stmt) {
    return evalString("(function() { try { " + stmt + "; return 'ok'; } catch (e) {\n"
                      "  return 'denied|' + e.message + '|' + e.code + '|' + e.path; } })()");
}

std::string refusal(const fs::path& p) {
    const std::string path = gatePath(p);
    return "denied|no files access: " + path + "|EACCES|" + path;
}

}  // namespace

void bromeshTestAccessGate() {
    std::cout << "[gate] Loaders and savers under the access gate..." << std::endl;

    const fs::path root = fs::temp_directory_path() / "bromesh_gate";
    const fs::path ok = root / "ok";
    const fs::path hidden = root / "hidden";
    std::error_code ec;
    fs::remove_all(root, ec);
    fs::create_directories(ok);
    fs::create_directories(hidden);
    const std::string obj = "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n";
    const std::string rig = R"({"name":"gate_rig","bones":[{"name":"root","parent":""}]})";
    for (const fs::path& dir : {ok, hidden}) {
        std::ofstream(dir / "tri.obj", std::ios::binary) << obj;
        std::ofstream(dir / "rig.json", std::ios::binary) << rig;
    }

    const std::string hiddenPrefix = gatePath(hidden);
    std::atomic<int> calls{0};
    ev::setAccessGate([&](const ev::AccessRequest& r) -> std::string {
        if (r.kind != ev::AccessKind::Read && r.kind != ev::AccessKind::Write) return "";
        ++calls;
        const std::string path(r.subject);
        return path.compare(0, hiddenPrefix.size(), hiddenPrefix) == 0 ? "no files access: " + path : "";
    });

    std::string r;
    // Readers.
    r = syncJs("Mesh.loadOBJ(" + jsStr(hidden / "tri.obj") + ")");
    GCHECK(r == refusal(hidden / "tri.obj"), r);
    r = syncJs("if (Mesh.loadOBJ(" + jsStr(ok / "tri.obj") + ").triangleCount !== 1) throw new Error('obj')");
    GCHECK(r == "ok", r);
    r = syncJs("Mesh.loadSTL(" + jsStr(hidden / "x.stl") + ")");
    GCHECK(r == refusal(hidden / "x.stl"), r);
    r = syncJs("Rig.specFromFile(" + jsStr(hidden / "rig.json") + ")");
    GCHECK(r == refusal(hidden / "rig.json"), r);
    r = syncJs("if (bro.rigging.specFromFile(" + jsStr(ok / "rig.json") + ").name !== 'gate_rig') throw new Error('rig')");
    GCHECK(r == "ok", r);
    // A relative path is asked about as the absolute path it opens.
    {
        const fs::path cwd = fs::current_path();
        fs::current_path(hidden);
        r = syncJs("Mesh.loadOBJ('tri.obj')");
        fs::current_path(cwd);
        GCHECK(r == refusal(hidden / "tri.obj"), r);
    }

    // Writers: refused, nothing is written.
    r = syncJs("Mesh.box(1, 1, 1).saveOBJ(" + jsStr(hidden / "out.obj") + ")");
    GCHECK(r == refusal(hidden / "out.obj"), r);
    GCHECK(!fs::exists(hidden / "out.obj"), "a refused saveOBJ wrote");
    r = syncJs("Mesh.box(1, 1, 1).saveSTL(" + jsStr(hidden / "out.stl") + ")");
    GCHECK(r == refusal(hidden / "out.stl"), r);
    GCHECK(!fs::exists(hidden / "out.stl"), "a refused saveSTL wrote");
    r = syncJs("Mesh.saveSplatPLY(" + jsStr(hidden / "s.ply") +
               ", {positions: new Float32Array(3), scales: new Float32Array(3), rotations: new Float32Array(4),"
               " opacities: new Float32Array(1), sh: new Float32Array(3), shDegree: 0})");
    GCHECK(r == refusal(hidden / "s.ply"), r);
    GCHECK(!fs::exists(hidden / "s.ply"), "a refused saveSplatPLY wrote");
    r = syncJs("if (!Mesh.box(1, 1, 1).saveOBJ(" + jsStr(ok / "out.obj") + ")) throw new Error('save')");
    GCHECK(r == "ok", r);
    GCHECK(fs::exists(ok / "out.obj"), "an allowed saveOBJ did not write");

    // With the gate removed, nothing asks and everything works as before.
    ev::setAccessGate(nullptr);
    const int before = calls.load();
    r = syncJs("if (Mesh.loadOBJ(" + jsStr(hidden / "tri.obj") + ").triangleCount !== 1) throw new Error('obj')");
    GCHECK(r == "ok", r);
    r = syncJs("if (!Mesh.box(1, 1, 1).saveOBJ(" + jsStr(hidden / "out.obj") + ")) throw new Error('save')");
    GCHECK(r == "ok", r);
    GCHECK(fs::exists(hidden / "out.obj"), "an ungated saveOBJ did not write");
    GCHECK(calls.load() == before, "a removed gate was still asked");

    fs::remove_all(root, ec);
    std::cout << "  access gate [PASS]" << std::endl;
}
