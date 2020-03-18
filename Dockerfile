# Galaxy for BRIDGE
# VERSION       0.2
FROM bgruening/galaxy-stable, 19.05.1
LABEL \
  description="Galaxy BRIDGE" \
  maintainer="chrisbarnettster@gmail.com, snptha002@myuct.ac.za"

ENV GALAXY_CONFIG_BRAND BRIDGE

# Install Tools & Data types
#ADD bridgetoolbox_tools.yml $GALAXY_ROOT/tools.yaml
#RUN install-tools $GALAXY_ROOT/tools.yaml && \
#    /tool_deps/_conda/bin/conda clean --tarballs --yes > /dev/null && \
#    rm /export/galaxy-central/ -rf && \
#    mkdir -p $GALAXY_HOME/workflows


# --- This section to be removed soon, as tools get added to the Galaxy Toolshed #
# Install manual tools 
WORKDIR /galaxy-central
COPY BRIDGE /galaxy-central/tools/bridge
COPY galaxy/config/tool_conf.xml.sample /galaxy-central/tools/bridge/tool_conf.xml
RUN sed -i '/<toolbox/r /galaxy-central/tools/bridge/tool_conf.xml' /galaxy-central/config/tool_conf.xml.sample
COPY galaxy/config/datatypes_conf.xml.sample /galaxy-central/config/
COPY galaxy/lib/binary.py  /galaxy-central/lib/galaxy/datatypes/
COPY galaxy/lib/molecules.py /galaxy-central/lib/galaxy/datatypes/
COPY galaxy/config/galaxy.ini.sample /galaxy-central/config/
COPY galaxy/config/dependency_resolvers_conf.xml /galaxy-central/config/
ADD  executables /galaxy-central/tools/bridge/md_tools/
# Changing the ownership
RUN chown -R $GALAXY_USER:$GALAXY_USER /galaxy-central/tools/bridge
RUN chown $GALAXY_USER:$GALAXY_USER /galaxy-central/lib/galaxy/datatypes/binary.py
RUN chown $GALAXY_USER:$GALAXY_USER /galaxy-central/lib/galaxy/datatypes/molecules.py
RUN chown $GALAXY_USER:$GALAXY_USER /galaxy-central/config/galaxy.ini.sample
RUN chown $GALAXY_USER:$GALAXY_USER /galaxy-central/config/dependency_resolvers_conf.xml
# ---

# Styling
ADD welcome.html /etc/galaxy/web/welcome.html
ADD welcome.html /galaxy-central/static/welcome.html
ADD welcome.html /galaxy-central/static/welcome.html.sample

# Expose port 80 (webserver), 21 (FTP server), 8800 (Proxy)
EXPOSE :80
EXPOSE :21
EXPOSE :8800

# Autostart script that is invoked during container start
CMD ["/usr/bin/startup"]

