export const isRpresent = (molFile) => {
  let splittedMolFile = molFile.split("\n")
  //Remove first 3 lines to avoid having R in the header of the file
  let shorterMolFile = splittedMolFile.slice(2, -1)
  for (let i = 0; i < shorterMolFile.length; ++i ) {
    if (shorterMolFile[i].includes("R")) {
      return true
    }
  }
  return false
}
