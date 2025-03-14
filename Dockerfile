FROM ruby:slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    nodejs \
    git && \
    rm -rf /var/lib/apt/lists/*

# Clean up
RUN apt-get clean && \
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*  /tmp/* /var/tmp/*

# set the working directory
WORKDIR /workspace

# copy the Gemfile to the image
COPY Gemfile ./

# install jekyll and dependencies
RUN gem install --no-document jekyll bundler && bundle install --no-cache