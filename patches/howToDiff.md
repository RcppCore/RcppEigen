
### Hint

```sh
diff --exclude=CMakeLists.txt -ruw eigen-3.3.9/Eigen/ inst/include/Eigen/ > patches/eigen-3.3.9.diff
diff --exclude=CMakeLists.txt -ruw eigen-3.3.9/unsupported/Eigen/ inst/include/unsupported/Eigen/ >> patches/eigen-3.3.9.diff
```
