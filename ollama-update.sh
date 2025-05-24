#!/bin/env bash
# Fetch installed version
echo "checking ollama version:"

current=$(ollama --version | awk '{print $4}')
echo " installed version: $current"
if [[ $current != v* ]]; then current="v$current"; fi
# Fetch latest release tag from GitHub
latest=$(curl -s https://api.github.com/repos/ollama/ollama/releases/latest \
         | grep '"tag_name":' \
         | head -n1 \
         | cut -d '"' -f4)
if [[ $latest != v* ]]; then latest="v$latest"; fi
echo " latest release: $latest"
# Compare and update if needed
if [[ "$current" != "$latest" ]]; then
  echo "Updating Ollama: $current -> $latest"
  curl -fsSL https://ollama.com/install.sh | sh      # installs latest[2]
else
  echo "Ollama is up-to-date ($current)"
fi
echo "Checking local models for updates.."
sleep 1
ollama list | awk 'NR>1 {print $1}' | \
 xargs -I {} sh -c 'echo "Updating model: {}"; ollama pull {}; echo "--"' \
 && echo "All models updated."
