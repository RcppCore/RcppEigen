
### Hint

```sh
diff --exclude=CMakeLists.txt -ruw eigen-3.4.0/Eigen/ inst/include/Eigen/ > patches/eigen-3.4.0.diff
diff --exclude=CMakeLists.txt -ruw eigen-3.4.0/unsupported/Eigen/ inst/include/unsupported/Eigen/ >> patches/eigen-3.4.0.diff
```
