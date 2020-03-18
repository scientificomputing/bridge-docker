
# Biomolecular Reaction & Interaction Dynamics Environment (BRIDGE)

[BRIDGE](https://github.com/scientificomputing/BRIDGE) is a flavour of the Galaxy platform to run and analyze molecular dynamics simulations.

## Get from DockerHub
```
docker pull scientificomputing/bridge
docker run -d -p 8080:80 scientificomputing/bridge
```

For more run options see [docker-galaxy-stable](https://github.com/bgruening/docker-galaxy-stable).

### Tagged versions

In case the latest version is not working, also see the tagged versions that are available.
For example try 19.05.1 which is pinned to Galaxy 19.05.1

```
docker pull scientificomputing/bridge:19.05.1
docker run -d -p 8080:80 scientificomputing/bridge:19.05.1
```

To find tags it is best to browse to the [BRIDGE DockerHub page](https://hub.docker.com/r/scientificomputing/bridge). Alternatively try the following on the commandline to list the first few tags:
```
curl "https://hub.docker.com/v2/repositories/scientificomputing/bridge/tags?page=1" |jq '."results"[]["name"]'
```


## Quick start with custom Galaxy Docker
- Build: `docker build -t bridge .`
- Run: `docker run -d -p 8080:80 --rm bridge` # the --rm flag will automatically remove the container when it exits
- Use: Open your web browser on http://localhost:8080

## Further information

See [BRIDGE](https://github.com/scientificomputing/BRIDGE)
