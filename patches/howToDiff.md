
### Hint

```sh
diff --exclude=CMakeLists.txt -ruw eigen-3.4.0/Eigen/ inst/include/Eigen/ > patches/eigen-3.4.0.diff
diff --exclude=CMakeLists.txt -ruw eigen-3.4.0/unsupported/Eigen/ inst/include/unsupported/Eigen/ >> patches/eigen-3.4.0.diff
```

or when using a git checkout of eigen (at the appropriate tag and branch)

```sh
diff --exclude=CMakeLists.txt -ruw ../eigen/Eigen/ inst/include/Eigen > patches/eigen-5.0.1.diff
diff --exclude=CMakeLists.txt -ruw ../eigen/unsupported/Eigen/ inst/include/unsupported/Eigen/ >> patches/eigen-5.0.1.diff
```
