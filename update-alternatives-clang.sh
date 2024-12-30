#!/bin/bash
## might be easier adding to PATH: /usr/lib/llvm-$ver/bin
# 241229
# add clang ${version} to Ubuntu

update_alternatives() {
    local version=${1}
    local priority=${2}
    local master=${3}
    local slaves=${4}
    local path=${5}
    local cmdln

    cmdln="--verbose --install ${path}${master} ${master} ${path}${master}-${version} ${priority}"
    for slave in ${slaves}; do
        cmdln="${cmdln} --slave ${path}${slave} ${slave} ${path}${slave}-${version}"
    done
    sudo update-alternatives ${cmdln}
}


if [[ ${#} -ne 2 ]]; then
    echo usage: "${0}" clang_version priority
    exit 1
fi

version=${1}
priority=${2}
path="/usr/bin/"

sudo apt update

# download and launch the setup script
curl -sSL https://apt.llvm.org/llvm.sh | sudo bash -s ${version}
# curl -sSL https://apt.llvm.org/llvm.sh | sudo bash -s ${version} all

# configure with update-alternatives
master="llvm-config"
# up to v17
slaves="llvm-addr2line llvm-ar llvm-as llvm-bcanalyzer llvm-bitcode-strip llvm-cat llvm-cfi-verify llvm-cov llvm-c-test llvm-cvtres llvm-cxxdump llvm-cxxfilt llvm-cxxmap llvm-debuginfo-analyzer llvm-debuginfod llvm-debuginfod-find llvm-diff llvm-dis llvm-dlltool llvm-dwarfdump llvm-dwarfutil llvm-dwp llvm-exegesis llvm-extract llvm-gsymutil llvm-ifs llvm-install-name-tool llvm-jitlink llvm-jitlink-executor llvm-lib llvm-libtool-darwin llvm-link llvm-lipo llvm-lto llvm-lto2 llvm-mc llvm-mca llvm-ml llvm-modextract llvm-mt llvm-nm llvm-objcopy llvm-objdump llvm-opt-report llvm-otool llvm-pdbutil llvm-PerfectShuffle llvm-profdata llvm-profgen llvm-ranlib llvm-rc llvm-readelf llvm-readobj llvm-reduce llvm-remark-size-diff llvm-remarkutil llvm-rtdyld llvm-sim llvm-size llvm-split llvm-stress llvm-strings llvm-strip llvm-symbolizer llvm-tapi-diff llvm-tblgen llvm-tli-checker llvm-undname llvm-windres llvm-xray"
# v18 and v19
slaves="llvm-addr2line llvm-ar llvm-as llvm-bcanalyzer llvm-bitcode-strip llvm-cat llvm-cfi-verify llvm-cov llvm-c-test llvm-cvtres llvm-cxxdump llvm-cxxfilt llvm-cxxmap llvm-debuginfo-analyzer llvm-debuginfod llvm-debuginfod-find llvm-diff llvm-dis llvm-dlltool llvm-dwarfdump llvm-dwarfutil llvm-dwp llvm-exegesis llvm-extract llvm-gsymutil llvm-ifs llvm-install-name-tool llvm-jitlink llvm-jitlink-executor llvm-lib llvm-libtool-darwin llvm-link llvm-lipo llvm-lto llvm-lto2 llvm-mc llvm-mca llvm-ml llvm-modextract llvm-mt llvm-nm llvm-objcopy llvm-objdump llvm-opt-report llvm-otool llvm-pdbutil llvm-PerfectShuffle llvm-profdata llvm-profgen llvm-ranlib llvm-rc llvm-readelf llvm-readobj llvm-readtapi llvm-reduce llvm-remarkutil llvm-rtdyld llvm-sim llvm-size llvm-split llvm-stress llvm-strings llvm-strip llvm-symbolizer llvm-tblgen llvm-tli-checker llvm-undname llvm-windres llvm-xray"
# v20
slaves="llvm-addr2line llvm-ar llvm-as llvm-bcanalyzer llvm-bitcode-strip llvm-cat llvm-cfi-verify llvm-cgdata llvm-cov llvm-c-test llvm-ctxprof-util llvm-cvtres llvm-cxxdump llvm-cxxfilt llvm-cxxmap llvm-debuginfo-analyzer llvm-debuginfod llvm-debuginfod-find llvm-diff llvm-dis llvm-dlltool llvm-dwarfdump llvm-dwarfutil llvm-dwp llvm-exegesis llvm-extract llvm-gsymutil llvm-ifs llvm-install-name-tool llvm-jitlink llvm-jitlink-executor llvm-lib llvm-libtool-darwin llvm-link llvm-lipo llvm-lto llvm-lto2 llvm-mc llvm-mca llvm-ml llvm-modextract llvm-mt llvm-nm llvm-objcopy llvm-objdump llvm-opt-report llvm-otool llvm-pdbutil llvm-PerfectShuffle llvm-profdata llvm-profgen llvm-ranlib llvm-rc llvm-readelf llvm-readobj llvm-readtapi llvm-reduce llvm-remarkutil llvm-rtdyld llvm-sim llvm-size llvm-split llvm-stress llvm-strings llvm-strip llvm-symbolizer llvm-tblgen llvm-tli-checker llvm-undname llvm-windres llvm-xray"

update_alternatives "${version}" "${priority}" "${master}" "${slaves}" "${path}"

master="clang"
# up to v17
slaves="asan_symbolize bugpoint clang++ clang-cpp clangd count dsymutil FileCheck ld64.lld ld.lld llc lld lldb lldb-argdumper lldb-instr lldb-server lldb-vscode lld-link lli lli-child-target not obj2yaml opt sanstats split-file UnicodeNameMappingGenerator verify-uselistorder wasm-ld yaml2obj yaml-bench"
# v18
slaves="asan_symbolize bugpoint clang++ clang-cpp clangd count dsymutil FileCheck ld64.lld ld.lld llc lld lldb lldb-argdumper lldb-dap lldb-instr lldb-server lld-link lli lli-child-target not obj2yaml opt sanstats split-file UnicodeNameMappingGenerator verify-uselistorder wasm-ld yaml2obj yaml-bench"
# v19 and v20
slaves="asan_symbolize bugpoint clang++ clang-cpp clangd count dsymutil FileCheck ld64.lld ld.lld llc lld lldb lldb-argdumper lldb-dap lldb-instr lldb-server lld-link lli lli-child-target not obj2yaml opt reduce-chunk-list sanstats split-file UnicodeNameMappingGenerator verify-uselistorder wasm-ld yaml2obj yaml-bench"

update_alternatives "${version}" "${priority}" "${master}" "${slaves}" "${path}"

exit


# to uninstall a Clang version
# LLVM_VERSION=18
# sudo apt purge -y clang-${LLVM_VERSION} lldb-${LLVM_VERSION} lld-${LLVM_VERSION} clangd-${LLVM_VERSION} && sudo apt autoremove -y
# sudo update-alternatives --config clang
# sudo update-alternatives --config llvm-config

# to generate the list of slaves on a clean Ubuntu vm
sudo apt update && sudo apt full-upgrade -y && sudo apt autoremove -y && sudo apt autoclean
ls /usr/bin > usr_bin_orig
version=20
curl -sSL https://apt.llvm.org/llvm.sh | sudo bash -s "$version"
# curl -sSL https://apt.llvm.org/llvm.sh | sudo bash -s "$version" all
ls /usr/bin > usr_bin_new
diff usr_bin_new usr_bin_orig | awk -v ver="$version" '/^< (llvm-).*-'"$version"'$/ && !/llvm-config-'"$version"'$/ {gsub("-" ver,""); printf "%s ", $2} END {print ""}' > slaves_llvm
diff usr_bin_new usr_bin_orig | awk -v ver="$version" '/^< / && !/llvm-/ && !/clang-'"$version"'$/ && /-'"$version"'$/ {gsub("^< ",""); gsub("-" ver "$",""); printf "%s ", $0} END {print ""}' > slaves_clang
# check if an application is already a master in:
# /var/lib/dpkg/alternatives/
