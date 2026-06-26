#!/usr/bin/env bash
echo "logical (nproc): $(nproc)"
phys=$(lscpu -p=Core 2>/dev/null | grep -v '^#' | sort -u | wc -l)
[ "$phys" -ge 1 ] 2>/dev/null || phys=$(nproc)
echo "physical cores: $phys"
echo "RAM GB: $(free -g | awk '/^Mem:/{print $2}')"
