pipeline{
    agent {
        // Equivalent to "docker build -f Dockerfile --build-arg version=1.0.2 ./build/
        dockerfile {
            filename 'Dockerfile'
            dir 'docker'
            label 'scCaTCH'
            additionalBuildArgs  '--build-arg version=0.0.1'
        }
    }
}