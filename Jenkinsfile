buildDockerImage([imageName: "sccatch", dockerContext: "docker", dockerFile: "docker/Dockerfile", pushRegistryNamespace: "obenauf", pushRegistry: "docker.artifactory.imp.ac.at", testCmd: null, pushBranches:["main"]])

//pipeline{
//    agent {
//        // Equivalent to "docker build -f Dockerfile --build-arg version=1.0.2 ./build/
//        dockerfile {
//            filename 'Dockerfile'
//            dir 'docker'
//            label 'scCaTCH'
//            additionalBuildArgs  '--build-arg version=0.0.1'
//        }
//    }
//}