#!/bin/bash
cd ~/NGS_data_test
git config --global user.email "follm@iarc.fr"
git config --global user.name "Circle CI_$CIRCLE_PROJECT_REPONAME_$CIRCLE_BRANCH"
git add .
git status
git commit --allow-empty -m "Results from $CIRCLE_PROJECT_REPONAME:$CIRCLE_BRANCH CircleCI tests build $CIRCLE_BUILD_NUM"
git push origin master

gem install github_changelog_generator
cd ~/needlestack
github_changelog_generator IARCbioinfo/needlestack --bug-labels bug,"minor bug" --no-pull-requests --future-release v1.1
git config --global user.email "follm@iarc.fr"
git config --global user.name "Circle CI_$CIRCLE_PROJECT_REPONAME_$CIRCLE_BRANCH"
git add .
git status
git commit -m "Added CHANGELOG"
git push origin $CIRCLE_BRANCH

curl -H "Content-Type: application/json" --data "{\"source_type\": \"Branch\", \"source_name\": \"$CIRCLE_BRANCH\"}" -X POST https://registry.hub.docker.com/u/iarcbioinfo/needlestack/trigger/fcd71f24-0f1e-40b8-b76a-09eb5a39e924/
