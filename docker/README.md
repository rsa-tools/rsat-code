
# Regulatory Sequence Analysis Tools Docker container

Please visit the documentation at 
[installing-RSAT](https://rsa-tools.github.io/installing-RSAT)
to learn how to pull and run the RSAT Docker container.


## Dockerfiles

Dockerfiles for the containers at [biocontainers/rsat](https://hub.docker.com/r/biocontainers/rsat) 
can be found at https://github.com/BioContainers/containers/tree/master/rsat

Those are based on this [Dockerfile](./Dockerfile), which we use for testing and buildind standalone containers.

## Building a new release of RSAT docker

1. Pull the latest version of RSAT code

```
cd rsat-code
git pull
```

2. Set a tag for this docker release

```
export RSAT_DOCKER_VERSION=`date '+%Y-%m-%d'`
git tag $RSAT_DOCKER_VERSION
```


3. Build the docker image

Beware: this step can take one or several hours. 

```
# assuming you were already in rsat-code, go to the docker sub-directory
cd docker

# Build a new image from the scratch
docker build --tag rsat:$RSAT_DOCKER_VERSION --tag rsat:latest --force-rm --compress --no-cache .
```

**Note:** the optiojn `--no-cache` ensures that the docker build incoroporate the recent changes in case a docker image was previously bult at the same place.

5. Submit the new image to

```
docker tag rsat:$RSAT_DOCKER_VERSION eeadcsiccompbio/rsat:$RSAT_DOCKER_VERSION
docker push eeadcsiccompbio/rsat:$RSAT_DOCKER_VERSION
```

**Note:** the push requires a login and password with write authorization on this repo. 

