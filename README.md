# vcf-patch

## Installation

```bash
cget install -f ./requirements.txt
mkdir build; cd build
cmake -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake -DCMAKE_BUILD_TYPE=Release ..
make
make install
```

## Example Usage
```
vcf-patch input.bcf patch_file.bcf chr6:251519-387461 > patched_output.bcf 
```
