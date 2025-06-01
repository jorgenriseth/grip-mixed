FROM ubuntu: 22.04

ENV PATH="/opt/greedy-1.3.0-alpha-Linux-gcc64/bin:$PATH"
RUN wget https://sourceforge.net/projects/greedy-reg/files/Nightly/greedy-nightly-Linux-gcc64.tar.gz/download -O /opt/greedy.tar.gz \
  && tar xvfz greedy.tar.gz \
  && rm greedy.tar.gz
