# function to perform supervised cell typing given profiles
# all that's done is computing logliks
# outputs: logliks, clusters, probs
# there should be a workhorse that's not exported, that gets called by the user-facing clustering and supervised functions