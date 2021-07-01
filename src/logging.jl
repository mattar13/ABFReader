#This is the bot that will be used to send the logging messages
const env_location = "C:\\Users\\mtarc\\OneDrive\\Documents\\JuliaCode\\.env"
dotenv(env_location) #First we can load the .env file

#this variable will be added
function BotNotify(text::String; add_date = true)
     if add_date
          println("[$(now())]: $text")
          sendMessage(text = "[$(now())]: $text")
     else
          println(text)
          sendMessage(text = text) 
     end
     return nothing
end

function wait_for_response()
     #We are going to run the bot until a message has been recieved
     run_bot() do msg
          response = msg.message.text
          println(lowercase(response))
          if lowercase(response) == "exit"
               println("Exiting loop")
               #We might need to do something weird like throw an error
               @assert 1 == 2 "InterruptException"
          else
               println(response)
          end
     end
end


