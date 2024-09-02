
# Regulatory Sequence Analysis Tools Docker container

Please visit the documentation at 
[installing-RSAT](https://rsa-tools.github.io/installing-RSAT)
to learn how to pull and run the **stable** RSAT Docker container.


## Dockerfiles

Dockerfiles for the **stable** containers at [biocontainers/rsat](https://hub.docker.com/r/biocontainers/rsat) 
can be found at https://github.com/BioContainers/containers/tree/master/rsat

Those are based on this [Dockerfile](./Dockerfile), which we use for testing and buildind standalone **development** containers.


## Building and uploading a new RSAT Docker development container

1. Get a git clone with the latest version of RSAT code

```
git clone https://github.com/rsa-tools/rsat-code.git
cd rsat-code
```

2. Set a tag for this docker release

```
export RSAT_DOCKER_VERSION=`date '+%Y-%m-%d'`
echo RSAT_DOCKER_VERSION=$RSAT_DOCKER_VERSION
git tag $RSAT_DOCKER_VERSION
git push origin --tags
```

3. Update rsat branch in Dockerfile

**Attention:** you also need to edit the file `rsat-code/docker/Dockerfile` in order to update the rsat branch. 
An example below, for the version 2024-08-28c:

```
RUN git clone https://github.com/rsa-tools/rsat-code.git  --branch 2024-08-28c --single-branch
```

4. Build the docker image

Beware: this step can take over one hour. 

```
# assuming you were already in rsat-code, go to the docker sub-directory
cd docker

# Build a new image from the scratch
docker build --tag rsat:$RSAT_DOCKER_VERSION --tag rsat:latest --force-rm --compress --no-cache .
```

**Note:** the option `--no-cache` ensures that the docker build incoroporate the recent changes in case a docker image was previously bult at the same place.

5. Submit the new image to

```
docker tag rsat:$RSAT_DOCKER_VERSION eeadcsiccompbio/rsat:$RSAT_DOCKER_VERSION
docker push eeadcsiccompbio/rsat:$RSAT_DOCKER_VERSION
```

**Note:** the push requires a login and password with write authorization on this repo. 

6. Check the availability

You can chek the availability of RSAT containers on dockerhub at 
<https://hub.docker.com/r/eeadcsiccompbio/rsat/tags>




