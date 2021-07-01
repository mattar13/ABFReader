#This is the bot that will be used to send the logging messages
const env_location = "C:\\Users\\mtarc\\OneDrive\\Documents\\JuliaCode\\.env"
dotenv(env_location) #First we can load the .env file

#this variable will be added
function BotNotify(text::String; add_date = true)
     if add_date
          println("$(now()) $text")
          sendMessage(text = "$(now()) $text")
     else
          println(text)
          sendMessage(text = text) 
     end
     return nothing
end