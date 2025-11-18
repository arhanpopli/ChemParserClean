import isEmpty from 'lodash/isEmpty'
import isNil from 'lodash/isNil'
import { copyToClipboard } from 'quasar'

export const copyData = (data) => {
  if (isEmpty(data)) {
    let error = "Cannot copy empty text area"
    return error
  }
  else {
    copyToClipboard(data)
      .then(() => {
      })
      .catch(() => {

      })
  }
}
