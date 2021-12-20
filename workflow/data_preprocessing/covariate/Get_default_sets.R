# Modified from: https://github.com/RTIInternational/biocloud_docker_tools/blob/master/peer/v1.3/run_peer.R

# Original Author: Francois Aguet

library(peer, quietly=TRUE)  # https://github.com/PMBio/peer

# Start analysis:
model <- PEER()


# Get the default values:
paste0("PriorAlphaA:   ",PEER_getPriorAlphaA(model))
paste0("PriorAlphaB:   ",PEER_getPriorAlphaB(model))
paste0("PriorEpsA:     ",PEER_getPriorEpsA(model))
paste0("PriorEpsB:     ",PEER_getPriorEpsB(model))
paste0("Tolerance:     ",PEER_getTolerance(model))
paste0("VarTolerance:  ",PEER_getVarTolerance(model))

