import { api } from 'boot/axios'

const useMakeRequest = async ( url, data ) => {
  let error = null
  let response_data = null
  let isLoading;
  url = api.defaults.baseURL + url

  console.log(url)

  try {

    const response = await api.post(url, data)
    response_data = response.data
    let isLoading = false
  }
  catch(error){
      error = "There was an error: " + error
      let isLoading = false
    }

  return { response_data, error, isLoading }
}

export { useMakeRequest }
