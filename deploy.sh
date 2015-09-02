#!/bin/bash
cd ~/NGS_data_test
git config --global user.email "follm@iarc.fr"
git config --global user.name "Circle CI_$CIRCLE_PROJECT_REPONAME_$CIRCLE_BRANCH"
git add .
git status
git commit --allow-empty -m "Results from $CIRCLE_PROJECT_REPONAME:$CIRCLE_BRANCH CircleCI tests build $CIRCLE_BUILD_NUM"
git push origin master
curl -H "Content-Type: application/json" --data '&#123;"docker_tag_name": "$CIRCLE_BRANCH"&#125;' -X POST https://registry.hub.docker.com/u/mfoll/robust-regression-caller/trigger/e85e32e7-508d-4971-b621-e1a0409f711a/
