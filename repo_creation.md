# how this repo as been created

### Clone old and new repo

```
git clone git@depot.biologie.ens.fr:rsat ens_rsat_code
git clone git@github.com:rsa-tools/rsat-code.git github_rsat_code
```

### Use filter-repo to filter the directory

```
cd ens_rsat_code
git remote remove origin  # no needed, but it is more secure (should not be done if we want to update this repo)
git filter-repo --invert-paths --path public_html/demo_files --path public_html/motif_databases --path public_html/sample_outputs --force --preserve-commit-hashes
```
we use **--invert-paths** to exclude path we have migrated in other repo

See https://github.com/newren/git-filter-repo


### Get the commit in the new repo


```
cd ..
cd github_rsat_code
```

- as we will have a lot of merged file (as repo as not been updated since a long time)
```
git config merge.renameLimit 9999
```

- do the merge
```
git remote add ens_updated  ../ens_rsat_code
git fetch ens_updated
git merge ens_updated/master -s recursive -X theirs
```


git push -u origin master
