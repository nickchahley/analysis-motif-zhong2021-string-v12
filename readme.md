Need to initialize submodule(s) after cloning: `git submodule update --init`

The pipeline, as is, is a series of .py files defined in `dvc.yaml`. 

Pulling data requires credentials to access our s3 buckets: `.dvc` files are pointers to data files store on a remote (defined in `.dvc/config`) which can be pulled using `dvc pull` 

Pipeline can be executed using `dvc repro`. Individual stages with `dvc repro <stage>`. Either of those accept a `--force` flag if they don't cooperate.
