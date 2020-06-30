## -----------------------------------------------------------------------------
library('TreeTools', quietly = TRUE, warn.conflicts = FALSE)
library('TreeDist')
treesMatchingSplit <- c(
  AB.CDEF = TreesMatchingSplit(2, 4),
  ABC.DEF = TreesMatchingSplit(3, 3)
)
treesMatchingSplit

proportionMatchingSplit <- treesMatchingSplit / NUnrooted(6)
proportionMatchingSplit

splitInformation <- -log2(proportionMatchingSplit)
splitInformation

treesMatchingBoth <- TreesConsistentWithTwoSplits(6, 2, 3)
combinedInformation <- -log2(treesMatchingBoth / NUnrooted(6))

sharedInformation <- sum(splitInformation) - combinedInformation
sharedInformation

# Or more concisely:
SplitSharedInformation(n = 6, 2, 3)

## ----mackay-8-1, echo=FALSE, fig.width=4, out.width='50%', fig.height=3, fig.align='center'----
library('TreeDist')
H <- function(inBracket) {
  expression(paste(italic('H'), plain('('), italic(inBracket), plain(')')))
}
oldPar <- par(mar = c(3.1, 0.1, 0, 0.1))
joint <- Entropy(c(1, 2, 1, 2) / 6)
plot(0, type = 'n', xlim = c(0, joint), ylim = c(5, 0), axes = FALSE)
axis(1, at = c(0, 0.5, 1, 1.5, round(joint, 2)))
mtext('Entropy / bits', 1, 2)
rect(joint - 1, 3.1, 1, 3.9, col = "#56B4E9")
rect(0, 0.1, joint - 1, 0.9, col = "#F0E442", border = NA)
rect(1, 1.1, joint, 1.9, col = "#F0E442", border = NA)
rect(joint - 1, 1.1, 1, 1.9, col = "#56B4E9", border = NA)
rect(joint - 1, 0.1, 1, 0.9, col = "#56B4E9", border = NA)
text(1, 3.5, pos=4,
     expression(paste(italic(I), plain('('), italic('A;B'), plain(')'))))



rect(0, 2.1, joint, 2.9)
text(joint / 2, 2.5, 
     expression(paste(italic('H'), plain('('), italic('A, B'), plain(')'))))

rect(0, 0.1, 1, 0.9)
text(0.5, 0.5, 
     expression(paste(italic('H'), plain('('), italic(A), plain(')'))))


rect(joint - 1, 1.1, joint - 0, 1.9)
text(joint - 0.5, 1.5, 
     expression(paste(italic('H'), plain('('), italic(B), plain(')'))))


rect(0, 4.1, joint - 1, 4.9, col = "#F0E442")
rect(1, 4.1, joint, 4.9,     col = "#F0E442")
text(joint / 2, 4.5, 
     expression(paste(italic('H'['D']), plain('('), italic('A, B'), plain(')'))))
par(oldPar)

