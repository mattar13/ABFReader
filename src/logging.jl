#This is the bot that will be used to send the logging messages
const env_location = "C:\\Users\\mtarc\\OneDrive\\Documents\\JuliaCode\\.env"
dotenv(env_location) #First we can load the .env file

#this variable will be added
function notify(text::String) 
     println(text)
     sendMessage(text = text)
end